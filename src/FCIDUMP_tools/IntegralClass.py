from typing import Union

class FCIDUMPReader:
    def __init__(self, filename: str):
        """Initialize the FCIDUMP reader with a filename."""
        self.filename = filename
        self.one_body_integrals: dict[tuple[int, int], float] = {}  # (i,j) -> value
        self.two_body_integrals: dict[tuple[int, int, int, int], float] = {}  # (i,j,k,l) -> value
        self.orbital_energies: dict[int, float] = {}    # i -> value
        self.core_energy = 0.0
        self.header_lines: list[str] = []  # Store header lines for later use
        self.read_file()
        # integral dicts of the natural (original) ordering
        self._one_body_integrals_natural = self.one_body_integrals.copy()
        self._two_body_integrals_natural = self.two_body_integrals.copy()
        self._orbital_energies_natural = self.orbital_energies.copy()

    def read_file(self):
        """Read the FCIDUMP file and process the integrals."""
        with open(self.filename, 'r') as f:
            for line in f:
                self.header_lines.append(line)
                l = line.replace(',', ' ').replace('=', ' ').split()
                if 'NORB' in line:
                    self.norb = int(l[l.index('NORB') + 1])
                if 'NELEC' in line:
                    self.nelec = int(l[l.index('NELEC') + 1])
                if 'MS2' in line:
                    self.ms2 = int(l[l.index('MS2') + 1])
                if ('/'  in line) or ('END' in line):
                    break

            # Process the integral lines
            for line in f:
                parts = line.strip().split()
                if len(parts) != 5:
                    raise ValueError(f"Invalid line format: {line}")

                # Parse the line
                value = float(parts[0])
                i, j, k, l = map(int, parts[1:])

                # Categorize and store the integral
                self._store_integral(value, i, j, k, l)

    def _store_integral(self, value: float, i: int, j: int, k: int, l: int) -> None:
        """Store integral in the appropriate category based on indices."""
        if i > 0 and j > 0 and k > 0 and l > 0:
            # Two-body integral
            # Store in canonical order to respect 8-fold symmetry
            indices = self._canonical_indices(i, j, k, l)
            self.two_body_integrals[indices] = value
        elif i > 0 and j > 0 and k == 0 and l == 0:
            # One-body integral
            # Store in canonical order (i≤j)
            if i <= j:
                self.one_body_integrals[(i, j)] = value
            else:
                self.one_body_integrals[(j, i)] = value
        elif i > 0 and j == 0 and k == 0 and l == 0:
            # Orbital energy
            self.orbital_energies[i] = value
        elif i == 0 and j == 0 and k == 0 and l == 0:
            # Core energy
            self.core_energy = value

    def _canonical_indices(self, i: int, j: int, k: int, l: int) -> tuple[int, int, int, int]:
        """
        Return canonical ordering of indices to respect 8-fold symmetry:
        ijkl = jilk = ijlk = ... = klij = ...
        
        The canonical ordering is chosen as the one where:
        1. max(i,j) is compared with max(k,l)
        2. The pair with larger maximum comes first
        3. Within each pair, the larger index comes first
        """
        # Sort pairs internally (larger first) - this is more efficient
        ij = (max(i, j), min(i, j))
        kl = (max(k, l), min(k, l))

        # Use tuple comparison which naturally does lexicographic ordering
        # If ij >= kl, put ij first; otherwise put kl first
        if ij >= kl:
            return (ij[0], ij[1], kl[0], kl[1])
        else:
            return (kl[0], kl[1], ij[0], ij[1])

    def get_integral(self, i: int, j: int, k: int, l: int) -> float:
        """
        Get the integral value for given indices.
        For two-body integrals, provide all four indices.
        For one-body integrals, provide i, j and set k=l=0.
        For orbital energies, provide i and set j=k=l=0.
        For core energy, set i=j=k=l=0.
        """
        if i > 0 and j > 0 and k > 0 and l > 0:
            # Two-body integral with 8-fold symmetry
            indices = self._canonical_indices(i, j, k, l)
            return self.two_body_integrals.get(indices, 0.0)
        elif i > 0 and j > 0 and k == 0 and l == 0:
            # One-body integral
            if i <= j:
                return self.one_body_integrals.get((i, j), 0.0)
            else:
                return self.one_body_integrals.get((j, i), 0.0)
        elif i > 0 and j == 0 and k == 0 and l == 0:
            # Orbital energy
            return self.orbital_energies.get(i, 0.0)
        elif i == 0 and j == 0 and k == 0 and l == 0:
            # Core energy
            return self.core_energy
        else:
            raise ValueError("Invalid indices for integral retrieval.")

    def get_integral_dict(self, norb: int) -> dict[tuple[int, int, int, int], float]:
        """
        Returns a dictionary of all integrals with non-redundant indices.

        Args:
            norb: Number of orbitals (maximum index value)

        Returns:
            dict: Dictionary with keys as tuples (i,j,k,l) and values as integral values.
                 Only includes canonically ordered indices to avoid redundancy due to 8-fold symmetry.
        """
        result = {}

        # Core energy
        result[(0, 0, 0, 0)] = self.core_energy

        # Orbital energies (diagonal one-body terms)
        for i in range(1, norb + 1):
            value = self.orbital_energies.get(i, 0.0)
            result[(i, 0, 0, 0)] = value

        # One-body integrals (i,j,0,0)
        for i in range(1, norb + 1):
            for j in range(i, norb + 1):  # i <= j to avoid redundancy
                key = (i, j)
                value = self.one_body_integrals.get(key, 0.0)
                result[(i, j, 0, 0)] = value

        # Two-body integrals (i,j,k,l)
        # Generate only canonically ordered indices to avoid redundancy
        for i in range(1, norb + 1):
            for j in range(1, i + 1):  # j <= i
                for k in range(i, norb + 1):  # k >= i
                    if k == i:
                        l_start = j  # When k=i, l starts from j to avoid redundancy
                    else:
                        l_start = 1
                    for l in range(l_start, k + 1):  # l <= k
                        indices = self._canonical_indices(i, j, k, l)
                        value = self.two_body_integrals.get(indices, 0.0)
                        result[indices] = value

        return result

    def permute_integrals(self, permutation: tuple[int, ...], t_passive: bool = True) -> None:
        """Permute all integral dictionaries according to the given permutation."""
        # permutation is applied to the natural ordering
        if self._one_body_integrals_natural:
            self.one_body_integrals = self.permute_integral_dict(  # type: ignore[assignment, arg-type]
                self._one_body_integrals_natural, permutation, t_passive)

        if self._two_body_integrals_natural:
            self.two_body_integrals = self.permute_integral_dict(  # type: ignore[assignment, arg-type]
                self._two_body_integrals_natural, permutation, t_passive)

        if self._orbital_energies_natural:
            self.orbital_energies = self.permute_integral_dict(  # type: ignore[assignment, arg-type]
                self._orbital_energies_natural, permutation, t_passive)

    def permute_integral_dict(
        self, 
        integral_dict: dict[Union[tuple[int, ...], int], float], 
        permutation: tuple[int, ...], 
        t_passive: bool = True
    ) -> dict[Union[tuple[int, ...], int], float]:
        """
        Apply orbital permutation to integrals and maintain canonical ordering.
        
        Args:
            integral_dict: Dictionary of integrals with tuple or int keys
            permutation: Tuple representing the orbital permutation
            t_passive: If True, permutation uses passive notation; if False, active notation
        
        Returns:
            dict: New dictionary with permuted canonical indices
        """
        if not t_passive:
            permutation = (0, ) + tuple([permutation.index(i + 1) + 1 
                                        for i in range(len(permutation))])
        else:
            permutation = (0, ) + permutation

        # Determine integral type by inspecting first key
        key = next(iter(integral_dict))
        
        if isinstance(key, int):
            # Orbital energy (single integer key)
            orb_energy_permuted_dict: dict[int, float] = {}
            for i, value in integral_dict.items():  # type: ignore[misc]
                new_i = permutation.index(i)  # type: ignore[arg-type]
                orb_energy_permuted_dict[new_i] = value
            return orb_energy_permuted_dict  # type: ignore[return-value]
        elif isinstance(key, tuple):  # type: ignore[misc]
            key_len = len(key)
            if key_len == 4:
                # Two-body integral
                two_body_permuted_dict: dict[tuple[int, int, int, int], float] = {}
                for (i, j, k, l), value in integral_dict.items():  # type: ignore[misc]
                    new_i = permutation.index(i)
                    new_j = permutation.index(j)
                    new_k = permutation.index(k)
                    new_l = permutation.index(l)

                    new_indices = self._canonical_indices(new_i, new_j, new_k, new_l)
                    two_body_permuted_dict[new_indices] = value
                return two_body_permuted_dict  # type: ignore[return-value]
            elif key_len == 2:
                # One-body integral
                one_body_permuted_dict: dict[tuple[int, int], float] = {}
                for (i, j), value in integral_dict.items():  # type: ignore[misc]
                    new_i = permutation.index(i)
                    new_j = permutation.index(j)
                    if new_i <= new_j:
                        one_body_permuted_dict[(new_i, new_j)] = value
                    else:
                        one_body_permuted_dict[(new_j, new_i)] = value
                return one_body_permuted_dict  # type: ignore[return-value]
            else:
                raise ValueError(f"Invalid tuple length in integral dictionary: {key_len}")
        else:
            raise ValueError(f"Invalid key type in integral dictionary: {type(key)}")

    def summarize(self):
        """Print a summary of the integrals read from the file."""
        print(f"FCIDUMP file: {self.filename}")
        print(f"Core energy: {self.core_energy}")
        print(f"Number of orbital energies: {len(self.orbital_energies)}")
        print(f"Number of one-body integrals: {len(self.one_body_integrals)}")
        print(f"Number of two-body integrals: {len(self.two_body_integrals)}")

    def dump_integrals(self, filename: str, permutation: tuple[int, ...], digits: int = 14) -> None:
        """
        Write permuted integrals to an FCIDUMP file.
        
        Args:
            filename: Output filename
            permutation: Orbital permutation to apply
            digits: Number of decimal digits for formatting (default: 14)
        """
        self.permute_integrals(permutation)

        if not self.header_lines:
            print("No header lines to write.")
            return

        # For nice output formatting
        def space(value: float, total_width: int = 22) -> str:
            return ' ' * (total_width - len(f"{value:.{digits}f}"))

        with open(filename, 'w') as f:
            f.writelines(self.header_lines)
            if self.two_body_integrals:
                for (i, j, k, l), value in self.two_body_integrals.items():
                    f.write(f"{space(value)}{value:.{digits}f}    {i}    {j}    {k}    {l}\n")

            if self.one_body_integrals:
                for (i, j), value in self.one_body_integrals.items():
                    f.write(f"{space(value)}{value:.{digits}f}    {i}    {j}    0    0\n")

            if self.orbital_energies:
                for i, value in self.orbital_energies.items():
                    f.write(f"{space(value)}{value:.{digits}f}    {i}    0    0    0\n")

            f.write(f"{space(self.core_energy)}{self.core_energy:.{digits}f}    0    0    0    0\n")


# Example usage
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python IntegralClass.py <fcidump_file>")
        sys.exit(1)

    fcidump_file = sys.argv[1]
    reader = FCIDUMPReader(fcidump_file)
    reader.summarize()

    # Example of accessing integrals
    print("\nSome example integral values:")
    print(f"Core energy: {reader.get_integral(0, 0, 0, 0)}")

    if reader.orbital_energies:
        first_orb = next(iter(reader.orbital_energies.keys()))
        print(f"Orbital energy for orbital {first_orb}: {reader.get_integral(first_orb, 0, 0, 0)}")

    if reader.one_body_integrals:
        first_one_body = next(iter(reader.one_body_integrals.keys()))
        i, j = first_one_body
        print(f"One-body integral ({i},{j}): {reader.get_integral(i, j, 0, 0)}")

    if reader.two_body_integrals:
        first_two_body = next(iter(reader.two_body_integrals.keys()))
        i, j, k, l = first_two_body
        i,j,k,l = 1,2,3,4
        print(f"Two-body integral ({i},{j},{k},{l}): {reader.get_integral(i, j, k, l)}")

    integral_dict = reader.get_integral_dict(8)
    print(len(integral_dict))
