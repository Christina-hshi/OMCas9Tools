# OMCas9Tools
Computational tools for OM with CRISPR Cas9 labeling. It also supports digestion enzymes (e.g., DLE), and a mixture of digestion enzymes and Cas9 probes.

Try the pre-built binary under `release/src/OMCas9Tools`.
It was built on a x86_64 machine. If it cannot be run on your machine, please compile  by yourself following our instructions below.

## Dependencies
- boost (v1.80.0)
  <br>You are recommended to install it with `conda`.

## Compile & Install
The program is implemented in C++. We use `CMake` to facilitate the compilation process. Please make sure you have `CMake` installed on your machine.

For compilation, simply do
```
./build.sh
```
or you can build step by step by yourself.
```
mkdir -p build
cd build
cmake ..
make OMCas9Tools
```
For installation, you may add the directory of compiled binary to the system environment variable `PATH`, or move it to anywhere the `PATH` knows. Or you can call the program with its full/absolute path without installation.

To add it to the system environment variable `PATH`, do (in bash shell)
```
export PATH=$PATH:<Directory to your OMCas9Tools>
```

## Usages
`OMCas9Tools` is a set of tools. Run `OMCas9Tools` directly for more instructions, which is shown below.
```
./OMCas9Tools <command> <options>
<command> available:
  insilicoDigest        In silico digestion of DNA sequences.
  find20merNGG          Find all 20merNGG in DNA.
  analyze               For analysis only.
#Specify <command> to see <options>.
```

### insilicoDigest
It supports both restriction enzymes (e.g., DLE) and Cas9 probes (for which the PAM sequence is required), and also a mixture of them. Specially for Cas9 probes, user can set the parameters such that the program will search for targets with the number of mismatches up to a user given value. For more options, please see blow.
```
./OMCas9Tools insilicoDigest <options>
Options:
  -h [ --help ]                        print help messages
  -e [ --enzymes ] arg                 file containing enzyme info with one
                                       enzyme per line in TSV format: <Name>
                                       <Seq> <Digest_site> <Mismatch_max>
                                       <Type(1:Restriction Enzyme
                                       2:CRISPR/Cas9)> <Channel>
  -f [ --seqfiles ] arg                file containing sequence files with one
                                       file per line.
  -m [ --mismatch ] arg (=0)           maximum number of mismatches allowed
                                       when searching for targets during in
                                       silico digestion. [Deprecated! please
                                       set for each enzyme in enzyme file.]
  -d [ --mergeDis ] arg (=1000)        merge nearby labels <= <mergeDis>
  -o [ --output ] arg (=digested.cmap) file to save the digested map in CMAP
                                       format
  -t [ --thread ] arg (=1)             number of threads to run program.
```
Please note the multi-threading feature is on our TODO list, and is not supported yet, therefore you can ignore the `-t` option for now.

### find20merNGG
It is for finding all 20merNGG in the DNA sequences. The unknown regions with base `n` or `N` will be skipped. You can run `./OMCas9Tools find20merNGG` to see more options.
```
./OMCas9Tools find20merNGG <options>
Options:
  -h [ --help ]                       print help messages.
  -f [ --seqfiles ] arg               file containing sequence files with one
                                      file per line.
  --skiploc                           not output the loci of 20merNGG (only
                                      counts) to save space.
  -o [ --output ] arg (=20merNGG.tsv) output file.
```

### analyze
This tool is currently for internal test only. Therefore please ignore it for now.

## License
This program is distributed under BSD 2-Clause License.

## Contact
For any questions, please create/post an issue in this Github repository.
