BioCantor paradigm
==================

The basis of the BioCantor paradigm is that objects are linked by flexible parent/child relationships; flexible
operations can be performed by jumping around the resulting hierarchy. In most cases, parent/child relationships are
used to establish the parent as the frame of reference for the location of the child, though parent/child relationships
are not required to define a coordinate system.

Three main object types populate this paradigm. :class:`Location` objects represent blocked features (which can
represent, say, a multi-exon transcript). :class:`Sequence` objects hold sequence data. :class:`Parent` objects define
parent/child relationships. The :class:`Parent` class is very flexible in order to accommodate different types of
relationships. For example, a :class:`Parent` object can optionally hold a pointer to a :class:`Sequence`, meaning that
sequence is the frame of reference for an object with that :class:`Parent`. A :class:`Parent` object can optionally
hold a :class:`Location`, meaning that is the location of the child relative to that parent. :class:`Parent` has
several optional parameters which enable different types of relationships and operations. See "Instantiating BioCantor
objects" for concrete examples.

Objects hold pointers to their own (optional) parent. :class:`Parent` objects do not hold pointers to children and can
be reused for multiple children. :class:`Location` objects, :class:`Sequence` objects, and :class:`Parent` objects can
all have parents. Multi-level hierarchies are established when :class:`Parent` objects have their own parents.

See the remaining documentation for detailed explanations of how to instantiate and use these classes.