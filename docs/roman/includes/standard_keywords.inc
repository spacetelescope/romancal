
Standard ASDF metadata
++++++++++++++++++++++
The following table lists the attributes that are *required* to be present
in all reference files.  The first column shows the attribute in the ASDF
reference file headers, which is the same as the name within the data
model meta tree (second column).  The second column gives the roman
data model name for each attribute, which is useful when using data
models in creating and populating a new reference file.

=========== ===============================
Attribute      Fully Qualified Path
=========== ===============================
author          model.meta.author
model_type      model.meta.model_type
date            model.meta.date
description     model.meta.description
instrument      model.meta.instrument.name
reftype         model.meta.reftype
telescope       model.meta.telescope
useafter        model.meta.useafter
=========== ===============================


**NOTE:** More information on standard required attributes can be found here:
:ref:`Standard ASDF metadata`
