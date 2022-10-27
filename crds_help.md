usage: /Users/mteodoro/.pyenv/versions/3.9.0/envs/my_dev_romancal_venv390/lib/python3.9/site-packages/crds/list.py
       [-h] [--references] [--mappings] [--cached-references]
       [--cached-mappings] [--cached-pickles] [--full-path]
       [--dataset-ids-for-instruments INSTRUMENTS [INSTRUMENTS ...]]
       [--dataset-headers IDS [IDS ...]] [--id-expansions-only]
       [--first-id-expansion-only] [--minimize-headers] [--json] [--config]
       [--status] [--cat [FILES ...]] [--keywords KEYWORDS [KEYWORDS ...]]
       [--add-filenames] [--no-arrays] [--tpns [FILES ...]]
       [--operational-context] [--remote-context PIPELINE]
       [--resolve-contexts] [--required-parkeys]
       [--file-properties [FILE_PROPERTIES ...]] [--contexts [CONTEXT ...] |
       --range MIN:MAX | --all | --last-n-contexts N | --up-to-context CONTEXT
       | --after-context CONTEXT] [--include-orphans] [-v]
       [--verbosity VERBOSITY] [--dump-cmdline] [-R] [-I] [-V] [-J] [-H]
       [--roman] [--stats] [--profile PROFILE] [--log-time] [--pdb]
       [--debug-traps]

crds.list prints out a variety of information about CRDS configuration, the
cache, reference and mapping files, default context names, and dataset headers
and ids used for CRDS reprocessing recommendations.


optional arguments:
  -h, --help            show this help message and exit
  --references          print names of reference files referred to by contexts
  --mappings            print names of mapping files referred to by contexts
  --cached-references   print the paths of all references in the local cache. very primitive.
  --cached-mappings     print the paths of all mappings in the local cache. very primitive.
  --cached-pickles      print the paths of all pickles in the local cache. very primitive.
  --full-path           print the full paths of listed files.
  --dataset-ids-for-instruments INSTRUMENTS [INSTRUMENTS ...]
                        print the dataset ids known to CRDS associated for the specified instruments.
  --dataset-headers IDS [IDS ...]
                        print matching parameters for the specified dataset ids.
  --id-expansions-only  print out only the <product>:<exposure> expansion associated with the specified --dataset-headers ids.
  --first-id-expansion-only
                        print out only the first exposure ID (header or expanded) associated with a particular product ID.
  --minimize-headers    print out only header parameters required by a particular CRDS context.
  --json                print out header parameters in JSON format suited for crds.bestrefs and grepping.
  --config              print CRDS configuration information.
  --status              print brief, basic, CRDS configuration information.
  --cat [FILES ...]     print the text of the specified mapping files or header and format info for references.
  --keywords KEYWORDS [KEYWORDS ...]
                        limited list of keywords to be catted from reference headers.
  --add-filenames       prefix each line of a cat'ed file with the filename.
  --no-arrays           Don't --cat array properties that are slow to compute. Use for large files.
  --tpns [FILES ...]    print the certify constraints (.tpn's) associated with the specified or implied files.
  --operational-context
                        print the name of the operational context on the central CRDS server.
  --remote-context PIPELINE
                        print the name of the context reported as in use by the specified pipeline.
  --resolve-contexts    print the literal names of the contexts defined by the command line context specifiers.
  --required-parkeys    print the names of the parkeys required to compute bestrefs for the specified mappings.
  --file-properties [FILE_PROPERTIES ...]
                        print the instrument, filekind, filename for each of the files specified.
  --contexts [CONTEXT ...]
                        Specify a list of CRDS mappings to operate on: .pmap, .imap, or .rmap or date-based specification
  --range MIN:MAX       Operate for pipeline context ids (.pmaps) between <MIN> and <MAX>.
  --all                 Operate with respect to all known CRDS contexts.
  --last-n-contexts N   Operate with respect to the last N contexts.
  --up-to-context CONTEXT
                        Operate on all contexts up to and including the specified context.
  --after-context CONTEXT
                        Operate on all contexts after and including the specified context.
  --include-orphans     Include reference files not mentioned by any contexts.
  -v, --verbose         Set log verbosity to True,  nominal debug level.
  --verbosity VERBOSITY
                        Set log verbosity to a specific level: 0..100.
  --dump-cmdline        Dump the command line parameters used to start the script to the log.
  -R, --readonly-cache  Don't modify the CRDS cache.  Not compatible with options which implicitly modify the cache.
  -I, --ignore-cache    Download required files even if they're already in the cache.
  -V, --version         Print the software version and exit.
  -J, --jwst            Force observatory to JWST for determining header conventions.
  -H, --hst             Force observatory to HST for determining header conventions.
  --roman               Force observatory to Roman for determining header conventions.
  --stats               Track and print timing statistics.
  --profile PROFILE     Output profile stats to the specified file.
  --log-time            Add date/time to log messages.
  --pdb                 Run under pdb.
  --debug-traps         Bypass exception error message traps and re-raise exception.

General categories of information driven by switches include:

0. Overall CRDS configuration
1. CRDS server file lists
2. CRDS cache file lists and paths
3. Cached file contents or headers
4. CRDS reprocessing dataset ids and parameters
5. Listing global default and installed pipeline contexts
6. Resolving context specifiers into literal context names

Many crds list services require setting CRDS_SERVER_URL to a valid CRDS
server to provide a source for the headers.

For HST::

% export CRDS_SERVER_URL=https://hst-crds.stsci.edu

or for JWST::

% export CRDS_SERVER_URL=https://jwst-crds.stsci.edu

0. Configuration information governing the behavior of CRDS for simple
configurations can be dumped::

crds list --status
CRDS Version = '7.0.7, bump-version, 7432326'
CRDS_MODE = 'auto'
CRDS_PATH = '/Users/jmiller/crds_cache_ops'
CRDS_SERVER_URL = 'https://jwst-crds.stsci.edu'
Effective Context = 'jwst_0204.pmap'
Last Synced = '2016-09-20 08:00:09.115330'
Python Executable = '/Users/jmiller/anaconda/bin/python'
Python Version = '3.5.2.final.0'
Readonly Cache = False

More comprehensive configuration information is also available for advanced
configurations::

% crds list --config
... lots of info ....

1. Files known by the CRDS server to belong to specified contexts can be listed
even if the files are not installed in a local CRDS Cache.

The --mappings command recursively evaluates and includes all the sub-mappings,
i.e. imaps and pmaps, of the specified contexts.

Contexts to list can be specified in a variety of ways:

-- To list the references contained by several contexts::

% crds list  --references --contexts hst_0001.pmap hst_0002.pmap ...
vb41935ij_bia.fits
vb41935kj_bia.fits
...

-- To list the references in a numerical range of contexts::

% crds list --references --range 1:2 --references
vb41935lj_bia.fits
vb41935oj_bia.fits
...

-- To list all mappings, even those not referenced by an imap or pmap::

% crds list --mappings --all
hst.pmap
hst_0001.pmap
hst_0002.pmap
hst_acs.imap
hst_acs_0001.imap
hst_acs_0002.imap
hst_acs_atodtab.rmap
...

--references, --mappings, or both can be listed.

2. Locally cached files (files already synced to your computer) can be listed::

% crds list --cached-mappings --full-path
...

% crds list --cached-references --full-path
...

In both cases adding --full-path prints the path of the file within the CRDS cache.

These are merely simple directory listings which ignore the context specifiers
and can be grep'ed for finer grained answers.

3. The contents of cached mappings or references (header only) can be printed to stdout like this::

% crds list --contexts jwst-fgs-linearity-edit jwst-nirspec-linearity-edit --cat --add-filename | grep parkey
CRDS - INFO - Symbolic context 'jwst-fgs-linearity-edit' resolves to 'jwst_fgs_linearity_0008.rmap'
CRDS - INFO - Symbolic context 'jwst-nirspec-linearity-edit' resolves to 'jwst_nirspec_linearity_0009.rmap'
/cache/path/mappings/jwst/jwst_fgs_linearity_0008.rmap:     'parkey' : (('META.INSTRUMENT.DETECTOR', 'META.SUBARRAY.NAME'), ('META.OBSERVATION.DATE', 'META.OBSERVATION.TIME')),
/cache/path/mappings/jwst/jwst_nirspec_linearity_0009.rmap:     'parkey' : (('META.INSTRUMENT.DETECTOR', 'META.SUBARRAY.NAME'), ('META.OBSERVATION.DATE', 'META.OBSERVATION.TIME')),

this prints the contents of the specified rmaps.

The -edit specifier above refers to mappings contained by the default starting point (.pmap) of future
server submissions.  It tracks on-going submission work that precedes the adoption of a new context
as the default in use by the pipeline.

crds.list --cat can be applied to references and prints out the reference metadata that CRDS views
abstractly as the file header.

References need to be catted explicitly by name,  but the list can come from the --references command
explained above::

% crds list --cat jwst_nirspec_dark_0036.fits
CRDS - INFO - Symbolic context 'jwst-operational' resolves to 'jwst_0167.pmap'
##########################################################################################
File:  '/grp/crds/jwst/references/jwst/jwst_nirspec_dark_0036.fits'
##########################################################################################
{'A1_COL_C': '8.9600000e+002',
'A1_CONF1': '2.1846000e+004',
...
}

4. Information about the dataset IDs and parameters used for CRDS reprocessing
and regressions can be printed or stored.

 Parameter set IDs can be listed for one or more instruments as follows::

 % crds list --dataset-ids-for-instruments wfc3...
 JCL403010:JCL403ECQ
 ... hundreds to hundreds of thousands of IDs as shown above ...

 IDs can also be captured to a file using UNIX I/O redirection::

 % crds list --dataset-ids-for-instruments wfc3   >wfc3.ids

 IDs for HST are of the form <product>:<exposure> where many exposures feed into
 the construction of one product and recalibrating any component exposure suggests
 recalibrating the combined product.

 CRDS stores dataset parameters for regression testing as a JSON dictionaries
 specifying one set of dataset parameters per line of the file::

 % crds list --dataset-headers @wfc3.ids --json > wfc3.headers.json

 NOTE:  while IDs can be specified directly on the command line,  CRDS has an
 @-notation that means "take IDs from this file".

 The JSON headers are suitable for running through crds.bestrefs to perform
 reprocessing checks or single context reference file coverage checks shown  here::

 % crds bestrefs --load-pickle wfc3.headers.json --dump-unique-errors --stats
 ...  errors related to looking up references for these parameter sets ...

 The script crds_dataset_capture combines the process of dumping all IDs for an
 instrument and dumping their corresponding dataset parameters.  IDs files and
 header files are placed in a dated regression capture directory::

 % crds_dataset_capture wfc3 acs ...
 ... downloads IDs and headers for WFC3, ACS to dated directory ...

 The default multi-line format for dataset parameters is more readable than the
 --json form::

 % crds list --dataset-headers jcl403010 --first-id --minimize-header
 CRDS - INFO - Symbolic context 'hst-operational' resolves to 'hst_0462.pmap'
 CRDS - INFO - Dataset pars for 'JCL403010:JCL403ECQ' with respect to 'hst_0462.pmap'
 {'APERTURE': 'WFC1',
  'ATODCORR': 'OMIT',
  'BIASCORR': 'COMPLETE',
  'CCDAMP': 'ABCD',
  'CCDCHIP': '-999.0',
  'CCDGAIN': '2.0',
  'CRCORR': 'OMIT',
  'DARKCORR': 'COMPLETE',
  'DATE-OBS': '2016-02-20',
  'DETECTOR': 'WFC',
  'DQICORR': 'COMPLETE',
  'DRIZCORR': 'COMPLETE',
  'FILTER1': 'CLEAR1L',
  'FILTER2': 'F814W',
  'FLASHCUR': 'LOW',
  'FLATCORR': 'COMPLETE',
  'FLSHCORR': 'OMIT',
  'FW1OFFST': '0.0',
  'FW2OFFST': '0.0',
  'FWSOFFST': '0.0',
  'GLINCORR': 'UNDEFINED',
  'INSTRUME': 'ACS',
  'LTV1': '0.0',
  'LTV2': '0.0',
  'NAXIS1': '4144.0',
  'NAXIS2': '4136.0',
  'OBSTYPE': 'IMAGING',
  'PCTECORR': 'UNDEFINED',
  'PHOTCORR': 'COMPLETE',
  'RPTCORR': 'UNDEFINED',
  'SHADCORR': 'OMIT',
  'SHUTRPOS': 'A',
  'TIME-OBS': '17:32:29.666665',
  'XCORNER': '0.0',
  'YCORNER': '0.0',
  'dataset_id': 'JCL403010:JCL403ECQ'}

 Sometimes it's desirable to know the individual exposures CRDS associates with a product id::

 % crds list --dataset-headers jcl403010 --id-expansions-only
 CRDS - INFO - Symbolic context 'hst-operational' resolves to 'hst_0462.pmap'
 JCL403010:JCL403ECQ
 JCL403010:JCL403EEQ
 JCL403010:JCL403EGQ
 JCL403010:JCL403EIQ
 JCL403010:JCL403EKQ
 JCL403010:JCL403EMQ
 JCL403010:JCL403EOQ
 JCL403010:JCL403EQQ
 JCL403010:JCL403ESQ
 JCL403010:JCL403EUQ

5. Information about the default context can be printed.  There are two variations and a subtle distinction::

% python m crds.list --operational-context
jwst_0204.pmap

lists the context which has been *commanded* as default on the CRDS server.

While::

% crds list --remote-context jwst-ops-pipeline
jwst_0101.pmap

lists the context which is *in actual use* in the associated archive pipeline as reported by
a cache sync echo.

During the interval between commanding a new default on the CRDS server and syncing the pipeline
CRDS cache,  the commanded and actual pipeline contexts can differ.

6. Resolving context specifiers

Some CRDS tools, including crds.list and crds.sync, support multiple
mechanisms for specifying context.  The --resolve-contexts command
interprets those specifiers into a non-recursive list of literal mapping
names and prints them out.  --resolve-contexts differs from --mappings
because it does not implicitly include all sub-mappings of the specified
contexts::

% crds list --resolve-contexts --all
jwst.pmap
jwst_0000.pmap
jwst_0001.pmap
jwst_0002.pmap
jwst_0003.pmap
...

% crds list --resolve-contexts --last 5
jwst_0205.pmap
jwst_0206.pmap
jwst_0207.pmap
jwst_0208.pmap
jwst_0209.pmap

% crds list --resolve-contexts  --contexts jwst-miri-dark-operational
jwst_miri_dark_0012.rmap

% crds list --resolve-contexts --contexts jwst-niriss-superbias-2016-01-01T00:00:00
jwst_niriss_superbias_0005.rmap
