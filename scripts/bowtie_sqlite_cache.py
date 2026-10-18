"""Bounded-memory alignment ingestion and exact coordinate-witness queries."""
import json
import sqlite3
from collections import Counter
from pathlib import Path

from bowtie_audit import parse_sam_line
from diffuprimer.io import sha256, utc_now


class AlignmentCache:
    def __init__(self, path):
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA journal_mode=WAL")
        self.connection.execute("PRAGMA synchronous=NORMAL")
        self.connection.execute("PRAGMA cache_size=-131072")
        self.connection.executescript("""
            CREATE TABLE IF NOT EXISTS oligos(id INTEGER PRIMARY KEY, sequence TEXT UNIQUE, rna_complete INT DEFAULT 0, genome_complete INT DEFAULT 0);
            CREATE TABLE IF NOT EXISTS subjects(id INTEGER PRIMARY KEY, name TEXT UNIQUE);
            CREATE TABLE IF NOT EXISTS hits(oligo INT, reference INT, subject INT, strand INT, start INT,
                PRIMARY KEY(oligo,reference,subject,strand,start)) WITHOUT ROWID;
            CREATE TABLE IF NOT EXISTS counts(oligo INT, reference INT, subject INT, strand INT, n INT,
                PRIMARY KEY(oligo,reference,subject,strand)) WITHOUT ROWID;
            CREATE TABLE IF NOT EXISTS imports(source TEXT PRIMARY KEY, path TEXT, reference INT, offset INT DEFAULT 0, complete INT DEFAULT 0);
            CREATE TABLE IF NOT EXISTS parts(left_oligo INT,right_oligo INT,reference INT,subject INT,start INT,end INT,
                PRIMARY KEY(left_oligo,right_oligo,reference,subject)) WITHOUT ROWID;
            CREATE TABLE IF NOT EXISTS parts_complete(left_oligo INT,right_oligo INT,reference INT,
                PRIMARY KEY(left_oligo,right_oligo,reference)) WITHOUT ROWID;
        """)
        self.connection.commit()
        self.subject_ids = dict((name, identifier) for identifier,name in self.connection.execute("SELECT id,name FROM subjects"))
        self.subject_names = {identifier:name for name,identifier in self.subject_ids.items()}

    def register(self, sequences):
        self.connection.executemany("INSERT OR IGNORE INTO oligos(sequence) VALUES(?)", [(s,) for s in sequences])
        self.connection.commit()
        return {sequence:(identifier,rna,genome) for identifier,sequence,rna,genome in self.connection.execute("SELECT id,sequence,rna_complete,genome_complete FROM oligos")}

    def ingest(self, path, queries, reference):
        source = str(reference) + ":" + sha256(path)
        row = self.connection.execute("SELECT offset,complete FROM imports WHERE source=?", (source,)).fetchone()
        if row and row[1]:
            return
        if not row:
            self.connection.execute("INSERT INTO imports(source,path,reference) VALUES(?,?,?)", (source,str(path),reference))
            self.connection.commit()
            row = (0,0)
        registered = self.register(queries.values())
        identifiers = {qid:registered[sequence][0] for qid,sequence in queries.items()}
        batch, counts, lines = [], Counter(), 0

        def flush(offset):
            self.connection.executemany("INSERT OR IGNORE INTO hits VALUES(?,?,?,?,?)", batch)
            self.connection.executemany("INSERT INTO counts VALUES(?,?,?,?,?) ON CONFLICT(oligo,reference,subject,strand) DO UPDATE SET n=n+excluded.n",
                [(*key,n) for key,n in counts.items()])
            self.connection.execute("UPDATE imports SET offset=? WHERE source=?", (offset,source))
            self.connection.commit()
            batch.clear()
            counts.clear()

        with Path(path).open("rb") as handle:
            handle.seek(row[0])
            while True:
                line = handle.readline()
                if not line:
                    break
                parsed = parse_sam_line(line.decode("ascii"), queries)
                if parsed is not None:
                    qid, hit = parsed
                    subject = hit["sseqid"]
                    if subject not in self.subject_ids:
                        cursor = self.connection.execute("INSERT INTO subjects(name) VALUES(?)", (subject,))
                        self.subject_ids[subject] = cursor.lastrowid
                        self.subject_names[cursor.lastrowid] = subject
                    item = (identifiers[qid],reference,self.subject_ids[subject],int(hit["strand"] == "-"),hit["start"])
                    batch.append(item)
                    counts[item[:4]] += 1
                lines += 1
                if lines % 50000 == 0:
                    flush(handle.tell())
                if lines % 5000000 == 0:
                    print(json.dumps({"ingestion": str(path), "reference": reference, "lines_since_resume": lines,
                                      "offset_bytes": handle.tell(), "at": utc_now()}), flush=True)
            flush(handle.tell())
        column = "rna_complete" if reference == 0 else "genome_complete"
        self.connection.executemany(f"UPDATE oligos SET {column}=1 WHERE id=?", [(identifier,) for identifier in identifiers.values()])
        self.connection.execute("UPDATE imports SET complete=1 WHERE source=?", (source,))
        self.connection.commit()

    def _part(self, left, right, reference, minimum=36, maximum=1000):
        if (minimum,maximum) != (36,1000):
            raise ValueError("Cache uses the frozen product-length rule")
        key = (left,right,reference)
        if not self.connection.execute("SELECT 1 FROM parts_complete WHERE left_oligo=? AND right_oligo=? AND reference=?",key).fetchone():
            left_length = self.connection.execute("SELECT length(sequence) FROM oligos WHERE id=?",(left,)).fetchone()[0]
            right_length = self.connection.execute("SELECT length(sequence) FROM oligos WHERE id=?",(right,)).fetchone()[0]
            plus = dict(self.connection.execute("SELECT subject,n FROM counts WHERE oligo=? AND reference=? AND strand=0",(left,reference)))
            minus = dict(self.connection.execute("SELECT subject,n FROM counts WHERE oligo=? AND reference=? AND strand=1",(right,reference)))
            for subject in sorted(plus.keys() & minus.keys()):
                witness = None
                if plus[subject] <= minus[subject]:
                    cursor = self.connection.execute("SELECT start FROM hits WHERE oligo=? AND reference=? AND subject=? AND strand=0",(left,reference,subject))
                    for (start,) in cursor:
                        lower,upper = start+max(left_length,minimum-right_length),start+maximum-right_length
                        found = self.connection.execute("SELECT start FROM hits WHERE oligo=? AND reference=? AND subject=? AND strand=1 AND start BETWEEN ? AND ? LIMIT 1",(right,reference,subject,lower,upper)).fetchone()
                        if found:
                            witness = (start,found[0]+right_length)
                            break
                else:
                    cursor = self.connection.execute("SELECT start FROM hits WHERE oligo=? AND reference=? AND subject=? AND strand=1",(right,reference,subject))
                    for (start,) in cursor:
                        lower,upper = start+right_length-maximum,min(start+right_length-minimum,start-left_length)
                        found = self.connection.execute("SELECT start FROM hits WHERE oligo=? AND reference=? AND subject=? AND strand=0 AND start BETWEEN ? AND ? LIMIT 1",(left,reference,subject,lower,upper)).fetchone()
                        if found:
                            witness = (found[0],start+right_length)
                            break
                if witness:
                    self.connection.execute("INSERT INTO parts VALUES(?,?,?,?,?,?)",(*key,subject,*witness))
            self.connection.execute("INSERT INTO parts_complete VALUES(?,?,?)",key)
            self.connection.commit()
        return [(subject,start,end) for subject,start,end in self.connection.execute("SELECT subject,start,end FROM parts WHERE left_oligo=? AND right_oligo=? AND reference=?",key)]

    def products(self, forward, reverse, reference):
        records = []
        for left_role,left in [("F",forward),("R",reverse)]:
            for right_role,right in [("F",forward),("R",reverse)]:
                for subject,start,end in self._part(left,right,reference):
                    records.append({"subject":self.subject_names[subject],"start":start,"end":end,
                        "amplicon_length":end-start,"left_primer":left_role,"right_primer":right_role})
        return records
