# MetaboAssayLibToMSPConversion
Conversion of OpenMS AssayGeneratorMetabo Assay Library to MSP.

Usage: parseMetaboAssayLibToMsp.py [OPTIONS]

Options:
  -in, --openmslib PATH           Input assay library from
                                  OpenMS::AssayGeneratorMetbo
                                  (.tsv,.traML,.pqp)
  -out, --msp PATH                Output spectral libray (.msp)
  --removedecoys / --no-removedecoys
                                  Removes the decoys from the assay library
                                  (default: True)
  --help                          Show this message and exit.


For question on how to use the utility please have a look at the 
conversion\_test, which can be used as utility test and an example.  
