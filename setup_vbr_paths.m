 % Download VBRc package
 % git clone https://github.com/vbr-calc/vbr.git
path_to_top_level_vbr='/path/to/vbr'; % point to vbr
addpath(path_to_top_level_vbr)

% Adds all relevant VBR paths to the matlab path
vbr_init

% Add functions directory to path
addpath('./functions/');