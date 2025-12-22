cwlVersion: v1.0
class: ExpressionTool

doc: |
  Create a subdirectory containing the specified files or directories
  This is used to organize workflow outputs into logical directories

requirements:
  - class: InlineJavascriptRequirement

inputs:
  dir_files:
    type:
      - type: array
        items: [File, Directory]
      - "null"
    doc: "Files or directories to include in the subdirectory"
  
  dir_name:
    type: string
    default: "output"
    doc: "Name of the subdirectory to create"

outputs:
  subdir:
    type: Directory
    doc: "Directory containing the input files"

expression: |
  ${
    // Handle different input types
    var items = [];
    
    if (inputs.dir_files == null) {
      items = [];
    } else {
      items = inputs.dir_files.filter(function(f) { return f != null; });
    }
    
    return {
      "subdir": {
        "class": "Directory",
        "basename": inputs.dir_name,
        "listing": items
      }
    };
  }
