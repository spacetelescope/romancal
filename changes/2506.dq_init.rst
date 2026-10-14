Skip adding the reference read and removing the data encoding offset in
``dq_init`` when the L1 metadata ``meta.instrument.l1_reference_read_added``
or ``meta.instrument.l1_data_encoding_offset_removed`` indicate these were
already done when the L1 file was made.
