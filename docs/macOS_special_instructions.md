# Special macOS instructions

A new conda install on macOS consistently causes the following error when first running samtools.

`dyld: Library not loaded: #rpath/libcrypto.1.0.0.dylib`

The fix is found [here](https://github.com/samtools/samtools/issues/974).

In summary... while in your non-base conda env run the following:

```bash
conda uninstall samtools
conda update --all
conda install samtools
conda install vsnp
```
