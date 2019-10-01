function [] = errorIfNotField( structure, field)
%CHECKISFIELD Function to check if field is part of structure, and give
%error message if not
if ~isfield(structure, field)
    error(['ERROR: ', field, ' is not set up.'])
end