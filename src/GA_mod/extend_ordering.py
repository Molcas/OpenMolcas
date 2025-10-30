# functions for extending a given ordering (chromosome).

def extend_ordering(ordering, on_site_permutation, num_prefix=0, num_suffix=0):
    """
    Extend a given ordering (chromosome) with a given prefix and suffix.

    Args:
        ordering: Tuple of integers representing the original permutation
        on_site_permutation: Tuple of integers representing the on-site pattern to repeat
        num_prefix: Number of prefix elements to add
        num_suffix: Number of suffix elements to add

    Returns:
        Tuple containing the extended permutation
    """

    n_onsite = len(on_site_permutation)

    # Create blocks for each position in ordering
    blocks = [
        [i * n_onsite + val for val in on_site_permutation]
        for i in range(len(ordering))
    ]

    # Arrange blocks according to ordering (subtract 1 to convert to 0-based index)
    tmp_ordering = tuple(val + num_prefix for i in ordering for val in blocks[i-1])
    prefix_tuple = tuple(range(1, num_prefix + 1))
    suffix_begin = num_prefix + len(tmp_ordering) + 1
    suffix_tuple = tuple(range(suffix_begin, suffix_begin + num_suffix))

    return prefix_tuple + tmp_ordering + suffix_tuple


def test_extend_ordering():
    """
    Example usage of extend_ordering function
    """
    # Example 1: Basic case
    ordering = (3, 1, 2)
    on_site = (1, 2, 4, 3)
    result = extend_ordering(ordering, on_site, 0, 0)
    print(f"Input ordering: {ordering}")
    print(f"On-site permutation: {on_site}")
    print(f"Extended ordering: {result}")
    # Output: (9,10,12,11,1,2,4,3,5,6,8,7)

    # Example 2: Different ordering
    ordering2 = (2, 1)
    on_site2 = (1, 3, 2)
    result2 = extend_ordering(ordering2, on_site2, 3, 2)
    print(f"\nInput ordering: {ordering2}")
    print(f"On-site permutation: {on_site2}")
    print(f"Extended ordering: {result2}")
    # Output: (1, 2, 3, 7, 9, 8, 4, 6, 5, 10, 11)

    # Example 3: Special on_site_permutation
    ordering3 = (5,1,2,4,3,7,6)
    on_site3 = (1,)
    suffix3 = 3
    result3 = extend_ordering(ordering3, on_site3, 0, 3)
    print(f"\nInput ordering: {ordering3}")
    print(f"On-site permutation: {on_site3}")
    print(f"The number of suffix: {suffix3}")
    print(f"Extended ordering: {result3}")

if __name__ == "__main__":
    test_extend_ordering()
