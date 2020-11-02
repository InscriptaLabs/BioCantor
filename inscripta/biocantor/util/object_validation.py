from inscripta.biocantor.exc import LocationException, LocationOverlapException


class ObjectValidation:
    @staticmethod
    def require_location_nonempty(location):
        if len(location) < 1:
            raise LocationException("Location must have positive length:\n{}".format(repr(location)))

    @staticmethod
    def require_location_has_parent(location):
        if not location.parent:
            raise LocationException("Location must have non-null parent attribute:\n{}".format(repr(location)))

    @staticmethod
    def require_location_has_parent_with_sequence(location):
        ObjectValidation.require_location_has_parent(location)
        if not location.parent.sequence:
            raise ValueError(
                "Parent of location must have non-null sequence attribute:\n{}".format(repr(location.parent))
            )

    @staticmethod
    def require_parent_has_location(parent):
        if not parent.location:
            raise ValueError("Parent must have non-null location attribute:\n{}".format(repr(parent)))

    @staticmethod
    def require_parent_has_parent(parent):
        if not parent.parent:
            raise ValueError("Parent must have non-null parent attribute:\n{}".format(repr(parent)))

    @staticmethod
    def require_parent_has_parent_with_location(parent):
        ObjectValidation.require_parent_has_parent(parent)
        if not parent.parent.location:
            raise ValueError(
                "Parent must have parent attribute with non-null location attribute:\n{}".format(repr(parent))
            )

    @staticmethod
    def require_parents_equal_except_location(parent1, parent2):
        if not parent1.equals_except_location(parent2):
            raise ValueError("Parents must be equal except location info:\n{}\n  !=\n{}".format(parent1, parent2))

    @staticmethod
    def require_locations_have_same_nonempty_parent(location1, location2):
        ObjectValidation.require_location_has_parent(location1)
        ObjectValidation.require_location_has_parent(location2)
        ObjectValidation.require_parents_equal_except_location(location1.parent, location2.parent)

    @staticmethod
    def require_locations_overlap(location1, location2, match_strand: bool = False):
        if not location1.has_overlap(location2, match_strand=match_strand):
            raise LocationOverlapException("Locations must overlap:\n{}\n{}".format(repr(location1), repr(location2)))

    @staticmethod
    def require_locations_do_not_overlap(location1, location2, match_strand: bool = False):
        if location1.has_overlap(location2, match_strand=match_strand):
            raise LocationOverlapException(
                "Locations must not overlap:\n{}\n{}".format(repr(location1), repr(location2))
            )

    @staticmethod
    def require_object_has_type(obj, required_type):
        if type(obj) is not required_type:
            raise TypeError("Object must have type {}".format(required_type))
