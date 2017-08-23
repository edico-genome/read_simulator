The read simulator supports is a modular framework that supports pipelines and modules to perform various simulation related tasks.

USAGE:
dSim.py -r run_config.json

AIM:
- Standardize our simulation environment, code re-use, minimize time to implement new (also unforeseen) requirements.

ARCHITECTURE:
- Framework supports simulation pipelines ( like germline, VLRD, ALT-Contig, Cancer )
- Simulation pipelines consist of different modules which are often re-used (e.g. bcftools/ Pirs).
- Common classes for
  Parsing and validating the config files
  Logging
  Interacting with the DB. 

PIPELINE EXAMPLE:
VLRD is a simulation pipeline that consists of 3 modules:
- custom variant simulator and truth VCF creator
- bcftools consensus mode to create Fasta
- Pirs

CREATING NEW MODULES
- A module (such as Pirs) is implemented as a Class that inherits from the ModuleBase class. The ModuleBase class specifies abstract methods (such as “validate config for my tool”,  “get_files_from_db”, “run”, and “update_db”). The abstract base methods force us to use consistent interfaces for all modules. The idea is that we can plug and play, and trivially replace e.g. my custom FASTA creator with BCFtools consensus tools or vice versa. A module is basically a wrapper around the tools we want to use. E.g. the Bcftools module invokes the Bcftools cmd line via the run method. The ModuleBase abstract class forces a consistent interface..

CREATING NEW MODULES
- A pipeline is also subclassed from the PipelinesBase class. This base class have fully defined methods that allows us to e.g. loop over its modules, validate the applicable configs setting, run all modules in pipeline, print results etc.

- As example, to create a new pipeline (VLRD) all we need to do is specify these 3 lines:
  // we inherit pipeline methods
  PipelineVLRD(extends PipelinesBase): 
 
 // specify an array of modules ( in order of desired execution ) that the pipeline should implement
 modules = [TruthVCFClass; FastaCreatorClass; Pirs]  
 
 // a pipeline Class is instantiated with it’s run settings
 PipelineVLRD constructor (run settings):
