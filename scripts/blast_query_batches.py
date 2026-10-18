"""Execute independent full-database queries without accumulating query batches."""
import json
import resource
import subprocess
from pathlib import Path

from diffuprimer.io import sha256,utc_now,write_json


def fasta(path):
    result={};name=None;parts=[]
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if name is not None:result[name]=''.join(parts)
            name=line[1:].split()[0];parts=[]
        else:parts.append(line)
    if name is not None:result[name]=''.join(parts)
    return result


def verify_tail(path):
    if not path.stat().st_size:return
    with path.open('rb') as handle:
        handle.seek(max(0,path.stat().st_size-4096));tail=handle.read()
    if not tail.endswith(b'\n') or len(tail.splitlines()[-1].split(b'\t'))!=14:
        raise ValueError('Incomplete BLAST alignment record')


def run_query_batches(command,run=subprocess.run,validate=True,batch_size=64):
    root=Path(__file__).resolve().parents[1]
    if validate:
        check=json.loads((root/'results/checks/blast_batching/check.json').read_text())
        if not check['exact_qualifying_hit_set_equality'] or check['script_sha256']!=sha256(Path(__file__)):
            raise ValueError('Query-batching equivalence check missing or stale')
    query=Path(command[command.index('-query')+1]);output=Path(command[command.index('-out')+1])
    database=Path(command[command.index('-db')+1])
    identity={'command':command,'query_sha256':sha256(query),'database_manifest_sha256':sha256(database.with_suffix('.manifest.json'))}
    directory=output.with_name(output.name+'.query_batches')
    directory.mkdir(parents=True,exist_ok=True)
    inputs=fasta(query);completed=[];missing=[]
    for index,(name,sequence) in enumerate(inputs.items()):
        mini_query=directory/f'{index:06d}.fasta';mini_output=directory/f'{index:06d}.tsv';marker=directory/f'{index:06d}.completed.json'
        mini_query.write_text(f'>{name}\n{sequence}\n')
        mini_command=command.copy()
        mini_command[mini_command.index('-query')+1]=str(mini_query)
        mini_command[mini_command.index('-out')+1]=str(mini_output)
        mini_identity={'command':mini_command,'queries_sha256':sha256(mini_query),'database_manifest_sha256':identity['database_manifest_sha256']}
        valid=False
        if marker.exists():
            prior=json.loads(marker.read_text())
            verify_tail(mini_output)
            valid=all(prior[k]==v for k,v in mini_identity.items()) and prior['output_sha256']==sha256(mini_output)
        if not valid:
            missing.append((index,name,sequence,mini_output,marker,mini_identity))
        completed.append(mini_output)
    groups=directory/'groups';groups.mkdir(exist_ok=True)

    def memory_limit():
        resource.setrlimit(resource.RLIMIT_AS,(32*1024**3,32*1024**3))

    def execute(batch):
        label=f'{batch[0][0]:06d}_{len(batch):03d}'
        grouped_query=groups/(label+'.fasta');grouped_output=groups/(label+'.tsv')
        grouped_query.write_text(''.join(f'>{r[1]}\n{r[2]}\n' for r in batch))
        actual=command.copy();actual[actual.index('-query')+1]=str(grouped_query);actual[actual.index('-out')+1]=str(grouped_output)
        with (groups/(label+'.stderr')).open('w') as handle:
            result=run(actual,check=False,stdout=subprocess.DEVNULL,stderr=handle,preexec_fn=memory_limit)
        if result.returncode:
            write_json(groups/(label+'.failed.json'),{'at':utc_now(),'command':actual,'exit_code':result.returncode,
                'address_space_limit_bytes':32*1024**3,'accepted_results':False})
            if len(batch)==1:raise subprocess.CalledProcessError(result.returncode,actual)
            middle=len(batch)//2
            execute(batch[:middle]);execute(batch[middle:]);return
        verify_tail(grouped_output)
        handles={r[1]:r[3].open('wb') for r in batch}
        try:
            with grouped_output.open('rb') as source:
                for line in source:
                    name=line.split(b'\t',1)[0].decode('ascii')
                    if name not in handles:raise ValueError('BLAST reported a query outside its input batch')
                    handles[name].write(line)
        finally:
            for handle in handles.values():handle.close()
        for index,name,sequence,path,marker,mini_identity in batch:
            verify_tail(path)
            write_json(marker,{**mini_identity,'actual_command':actual,'actual_queries_sha256':sha256(grouped_query),
                'output_sha256':sha256(path),'completed_at':utc_now(),
                'completion_evidence':'Checked independent-query batch exit 0; records separated by original query identifier'})
        print(json.dumps({'new_query_batch_complete':len(batch),'starting_query':batch[0][0],
            'total_queries':len(inputs),'reference':str(database),'at':utc_now()}),flush=True)

    for offset in range(0,len(missing),batch_size):execute(missing[offset:offset+batch_size])
    temporary=output.with_suffix(output.suffix+'.merging')
    with temporary.open('wb') as handle:
        for path in completed:
            with path.open('rb') as source:
                while block:=source.read(1024*1024):handle.write(block)
    temporary.replace(output)
    verify_tail(output)
    write_json(directory/'execution.json',{**identity,'completed_at':utc_now(),'queries':len(inputs),'output_sha256':sha256(output),
        'policy':'Each original query retains its identifier, sequence, full reference, scoring, E threshold, HSP/target ceilings and thread count; merge independent query outputs',
        'batch_size':batch_size,'address_space_limit_bytes':32*1024**3,
        'resource_failure_policy':'Unsuccessful batches contribute no accepted results; recursively split and repeat every query until a checked complete search succeeds'})
    return subprocess.CompletedProcess(command,0)
