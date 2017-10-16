# FragGeneScan++
An improved version of [FragGeneScan-Plus][fgsp].

## Installing
To build/install FragGeneScan++ from source, you will need [Meson][meson]. run the following commands from the source
repository:

```sh
meson build
ninja -C build
sudo ninja -C build install
```

## Running
An example usage of FragGeneScan++ looks like the following:

```
FGS++ -s example/NC_000913-454.fna -o output -w 0 -t 454_5 -p 16
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
