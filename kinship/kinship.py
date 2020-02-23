import numpy as np
import json


def load_population(directory='population/dominican.json'):
    # Open json file with frequencies
    pop_data = json.loads(open(directory, 'r').read())

    # Change key name for str(float) for downstream analyses
    for marker, freqs in pop_data.items():
        pop_data[marker] = {float(k): v for k, v in dict(freqs).items()}

    return pop_data


class Profile:
    """ STR profile from a unique individual """

    def __init__(self, profile: dict):
        """ Instance of Profile Class

        Args:
            alleles (dict): a dictionary with marker: alleles (list) value pairs.
            Example profile: {'AMEL': ['X', 'Y'], 'CSF1PO': [13, 14]}

        Returns:
            Profile object. The return code::

                Profile(markers=['AMEL', 'CSF1PO'],
                        alleles=[('X', 'Y'), (13, 14)])
        """
        self.profile = profile
        self.markers = self._parse_markers()
        self.alleles = self._parse_alleles()

    def _parse_markers(self):
        """ Returns a list of markers ordered lexicographically. """
        return sorted(self.profile.keys())

    def _parse_alleles(self):
        """ Returns a list of tuples of alleles ordered by marker name. """
        return np.array([tuple(self.profile[marker]) for marker in self.markers])

    def __repr__(self):
        return f'Profile(markers={self.markers}, alleles={self.alleles})'


class Duo:
    """ A paternity test for a Duo """

    def __init__(self, parent, child, population=load_population()):

        try:
            if sorted(parent.keys()) != sorted(child.keys()):
                raise ValueError()

        except ValueError:
            msg = """To calculate the LR of a hypothesis being correct
                both profiles must have the same markers."""
            raise ValueError(msg)

        self.child = self._load_profile(child)
        self.parent = self._load_profile(parent)
        self.markers = self._get_markers()
        self.population = population

    def _load_profile(self, data):
        try:
            if not isinstance(data, dict):
                msg = """
                Profile data must be a dictionary of type:
                {'marker': [allele_1, allele_2], ...} """
                raise TypeError(msg)

            for k, v in data.items():
                if len(v) != 2:
                    raise ValueError(f'Marker {k} must have two alleles')

        except:
            raise ValueError(
                'Something went wrong, check previous error messages.')

        if 'AMEL' in data.keys():
            del data['AMEL']

        return Profile(data)

    def _get_markers(self):
        try:
            if self.parent.markers != self.child.markers:
                msg = """Profiles introduced must have the same markers."""
                raise ValueError(msg)
        except:
            raise ValueError()

        return self.parent.markers

    @property
    def inconsistent_markers(self):
        """ Check if child contains at least one allele from its progenitor """

        _alleles = zip(self.parent.alleles, self.child.alleles)
        _markers, inconsistent_markers = self.parent.markers, []

        for marker, (p_alleles, c_alleles) in zip(_markers, _alleles):
            _check = [allele for allele in c_alleles if allele not in p_alleles]

            if len(_check) != 1:
                inconsistent_markers.append(marker)

        return inconsistent_markers if any(inconsistent_markers) else 0

    @property
    def paternity_index(self, calc=1):
        """ Calculate paternity index for Duo """

        if self.inconsistent_markers:
            return f'There are {len(self.inconsistent_markers)} inconsistent markers.'

        def is_heterozygous(alleles):
            return True if (len(np.unique(alleles)) == 2) else False

        for m, p, c in zip(self.markers, self.parent.alleles, self.child.alleles):
            pop = self.population[m]

            if is_heterozygous(p) and is_heterozygous(c):

                if len(np.intersect1d(p, c)) == 2:
                    a, b = pop[p[0]], pop[p[1]]
                    calc *= (a + b) / (4 * a * b)

                else:
                    shared_allele = np.intersect1d(p, c)[0]
                    calc *= 0.25 / pop[shared_allele]

            if is_heterozygous(p) and not is_heterozygous(c):
                calc *= 0.5 / pop[p[0]]

            if not is_heterozygous(p) and is_heterozygous(c):
                calc *= 0.5 / pop[p[0]]

            if not is_heterozygous(p) and not is_heterozygous(c):
                calc *= 1 / pop[pop[0]]

        return calc

    def __repr__(self):
        return f'Duo(parent={self.parent}, child={self.child})'
