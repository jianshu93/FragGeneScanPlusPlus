# FragGeneScan++

A modified fork of [FragGeneScan-Plus][fgsp].

**This tool is no longer maintained**. It was replaced by [FragGeneScanRs](https://github.com/unipept/FragGeneScanRs), a drop-in replacement for FragGeneScan, implemented in rust. FragGeneScanRs is both faster and contains a number of bug fixes that are present in both FragGeneScan and FragGeneScanPlusPlus.


## Installation

The installation instruction from [FGS+][fgsp install] remain valid:

1. Install `make` and `gcc`.

2. Clone this repository and go to the repository root.

   ```sh
   git clone https://github.com/unipept/FragGeneScanPlusPlus.git
   cd FragGeneScanPlusPlus
   ```

3. Run `make`.


## Running

An example usage of FragGeneScan++ looks like the following:

```
./FGSpp -s example/NC_000913-454.fna -o output -w 0 -t 454_5 -p 16
```

If you're running `FGSpp` from outside the repo, supply the training
directory location with `-r`:

```
FGSpp -s <this-repo>/example/NC_000913-454.fna -o output -w 0 -r <this-repo>/train -t 454_5 -p 16
```

For more info or a list of all options, run FragGeneScan++ without arguments.

### Output files

Upon completion, FragGeneScan++ can generate 3 files:

* The `[output_file].faa` file lists amino acid sequences in FASTA format
corresponding to the putative genes. The `[output_file].faa` file is always
generated.
* The `[output_file].ffn` file lists the nucleotide sequences in FASTA format
corresponding to the putative genes. This file is not automatically generated:
to obtain this file, set the `-d` flag to 1.
* The `[output_file].out` file lists the coordinates of putative genes. This
file consists of 5 columns ( _start position_, _end position_, _strand_,
_frame_, and _score_). This file is not automatically generated: to obtain this
file, set the `-e` flag to 1.


## License

FragGeneScan++ is released under under the terms of the GNU General Public
License, version 3 (or any later version), as published by the Free Software
Foundation.

Please see the [LICENSE](LICENSE) file for further information.


[fgsp]: https://github.com/hallamlab/FragGeneScanPlus/
[meson]: https://mesonbuild.com/
[fgsp install]: https://github.com/hallamlab/FragGeneScanPlus/wiki#setup-and-dependencies-1
