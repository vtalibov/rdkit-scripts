# pySMCM

[RDKit](http://www.rdkit.org/)-based Python script that calculates Synthethic Molecular Complexity Metrics (SMCM) descriptor, described by Allu&Oprea, *J. Chem. Inf. Model., 2005, 45, 1237-1243*.

### Example

```python
>>> from pySMCM import smcm_score
>>> smcm_score('CC(C)[C@@H](C(=O)N1CC2(CC2)C[C@H]1C3=NC=C(N3)C4=CC5=C(C=C4)C6=C(C5(F)F)C=C(C=C6)C7=CC8=C(C=C7)N=C(N8)[C@@H]9[C@H]1CC[C@H](C1)N9C(=O)[C@H](C(C)C)NC(=O)OC)NC(=O)OC')
130.51100000000002
```

### References

Electronegativity vaulues and SMCM score values for various bonds are taken form the original publication and reffered to Ivanciuc *et al.*, *J. Chem. Inf. Comput. Sci., 2001, 41, 269-272*.

SMARTS entries for scaffolds that improve chemical accessibility of molecules are taken from Chembank, 2005.
