# pickleprot_alpha

# Extra arguments in quick_load
Added an option to provide extra arguments into wildcard parsing. Extra arguments are prefaced with { (left curled
bracket) and separated by , (commas). For example, {2,3 in an extra argument would mean "use the first generator
to replace the first two wildcards with identical values, then use the second generator to replace three wildcards
with identical values". This solves issues when a filename in a dictionary provided by wildcards matches the name of the
folder, but is variable across different folders ( /path/model_2/structure_2.pdb ; /path/model_4/structure_4.pdb } which,
if handled by standard (in this case two) generators, would confuse the model loader by providing several non-existing 
paths.

EXAMPLES:

dict = { '/path/***/folder/***{2' : ( (1,2,3), 'some name' ) }

Will evaluate to:

/path/1/folder/1
/path/2/folder/2
/path/3/folder/3

dict = { 'path/***/folder/***/another/***/folder/***/***{2,3' : ( (1,2), ('a','b'), ('some other name') }

Will evaluate to:
path/1/folder/1/another/a/folder/a/a
path/1/folder/1/another/a/folder/a/a
path/2/folder/2/another/b/folder/b/b
path/2/folder/2/another/b/folder/b/b

several_new_models = proteinModel.quick_load(dict, some_handles_of_names)
