# Contribution Guidelines

## Tool Preparation for Test Toolshed

### Working Wrapper

The tool's wrapper should be valid and contain a basic set of parameters to get the tool running. In general, you should abide by the [IUC Best Practices](https://galaxy-iuc-standards.readthedocs.io). Wherever possible, use consistent terminology between tools (_e.g._ for parameter labels). 

### Auto-installation

The tool should be installable automatically. This implies that any dependencies should be wrapped into Galaxy packages, unless they already exist and are maintained by the devteam or IUC.

### Documentation

#### Citations

Appropriate citations should be included in the `citations.xml` file and imported into the wrapper (see the template tool for an example). This should include the latest reference for the Galaxy project (_e.g._ GCC presentation, BioRxiv preprint, or final publication) and references for any software tools used by the Galaxy tool.

#### README

A short description of the tool should be included in the wrapper's README section. This can be quoted from the underlying tool's website; make sure it is referenced properly. Whether or not you quote the description, you should include a link to the tool website for the user's convenience. 

### Toolshed Items

#### Package Synopsis (Tool Dependency Repository)

"Contains a tool dependency definition that downloads and compiles $tool version $tool_version"

$tool = RADIA, EXPANDS, Strelka, ProDuSe (proper caps)
$tool_version = 0.1.2, 1.0.14, etc.

#### Tool Synopsis

A single sentence explaining the purpose of the tool.

e.g. Identifies SNVs from tumour/normal BAMs 

#### Detailed Description

Not necessary, include at your leisure if you think text isn't redundant. 

## Tool Polishing for Main Toolshed

TBD.

