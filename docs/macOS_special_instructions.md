# Special macOS instructions

A new conda install on macOS consistently causes the following error when first running samtools.<br>
`dyld: Library not loaded: #rpath/libcrypto.1.0.0.dylib`

The fix is found here.<br>
https://github.com/samtools/samtools/issues/974

In summary... while in your non-base conda env run the following:<br>
`conda uninstall samtools`<br>
`conda update --all`<br>
`conda install samtools`<br>
`conda install vsnp`