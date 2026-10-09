# Notice and disclaimer

## Provenance

SIMPROT was developed by Paulo A. S. Nuin as a postdoctoral fellow in the laboratory of Elisabeth R. M. Tillier at the University Health Network (UHN), Toronto. The salary was paid from research grants. It is described in:

> Pang A, Smith AD, Nuin PAS, Tillier ERM (2005). SIMPROT: using an empirically determined indel distribution in simulations of protein evolution. *BMC Bioinformatics* 6:236. doi:10.1186/1471-2105-6-236

The original source (`legacy/simprot.cpp`) carries this notice, which is kept unchanged:

> Simprot (c) Copyright 2005-2012 by the University Health Network written by Elisabeth Tillier.
> Permission is granted to copy and use this program provided no fee is charged for it
> and provided that this copyright notice is not removed.

`legacy/random.c` contains code from `tools.c` in Ziheng Yang's PAML package, and it remains under PAML's terms. The Dayhoff and JTT data embedded in `simprot-cpp/tools/make_eigen.py` come from PAML's `dat/` files and from the publications cited there.

## License

The developer releases SIMPROT, including the C++20 reimplementation, tools and tests in `simprot-cpp/`, under the GNU General Public License, version 3 or (at your option) any later version (`LICENSE`), **to the extent that he holds rights in it**.

The ownership of the original code written at UHN has not been formally established. If UHN or another party holds rights in that code, those rights are unaffected by this release: the UHN notice above continues to apply to the original code, and, in particular, no fee may be charged for it. Anyone who holds such rights and objects to this release is asked to open an issue at https://github.com/nuin/simprot so that it can be resolved.

## No warranty

This software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and non-infringement. In no event shall the authors, the University Health Network or other copyright holders be liable for any claim, damages or other liability arising from the software or its use. Simulated data are not a substitute for empirical data, and results should be validated for any scientific or clinical purpose.
