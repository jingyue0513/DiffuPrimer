import torch

from diffuprimer.model import ConditionalDenoiser, corrupt, diffusion_loss, reveal_probability


def test_absorbing_terminal_and_fixed_region():
    labels = torch.tensor([[0, 1, 1, 0, 2, 2]])
    fixed = torch.tensor([[False, True, True, False, False, False]])
    noisy, mask = corrupt(labels, torch.tensor([16]), 16, fixed=fixed)
    assert noisy.tolist() == [[3, 1, 1, 3, 3, 3]]
    assert reveal_probability(1) == 1
    assert reveal_probability(16) == 1 / 16


def test_forward_loss_finite_and_backpropagates():
    model = ConditionalDenoiser(width=16, layers=2)
    dna = torch.randint(0, 5, (2, 64))
    clean = torch.zeros((2, 64), dtype=torch.long)
    clean[:, 3:23] = 1
    clean[:, 35:55] = 2
    valid = torch.ones_like(clean, dtype=torch.bool)
    t = torch.tensor([8, 16])
    noisy, masked = corrupt(clean, t, 16, valid)
    logits = model(dna, noisy, t / 16, torch.ones(2, 9), valid)
    loss = diffusion_loss(logits, clean, masked, t, 16, valid)
    assert torch.isfinite(loss)
    loss.backward()
    assert model.head[-1].weight.grad.abs().sum() > 0


def test_padding_does_not_contribute_to_loss():
    logits = torch.zeros(1, 4, 3)
    clean = torch.tensor([[1, 2, 0, 0]])
    valid = torch.tensor([[True, True, False, False]])
    masked = torch.ones_like(valid)
    loss = diffusion_loss(logits, clean, masked, torch.tensor([1]), 1, valid)
    assert abs(loss.item() - torch.log(torch.tensor(3.)).item()) < 1e-6
