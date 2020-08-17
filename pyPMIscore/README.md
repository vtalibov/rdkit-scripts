# pyPMIscore

[RDKit](http://www.rdkit.org/)-based Python script to evaluate 3D-alikeness of a small molecule.

The script calculates an average sum of normalized principal moments of intertia (I2/I1 + I3/I1) for a set of low-energy conformers. PMIscore value > 1.2 indicates some 3-dimentionality, > 1.6 significant 3-dimentionality.

### Example

```bash
$ python pyPMIscore.py 'O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl'
1.378836204079491
```

### References

[Prosser, Stokes and Cohen, ACS Med. Chem. Lett. 2020, 11, 1292-1298](https://pubs.acs.org/doi/10.1021/acsmedchemlett.0c00121)

[Sauer and Schwarz. J. Chem. Inf. Comput. Sci. 2003, 43, 3, 987-1003](https://pubs.acs.org/doi/10.1021/ci025599w)
