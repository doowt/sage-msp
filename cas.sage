from sage.all import *

def i_less_than_j(M, i, j):
    '''Compares two rows according to lexicographic ordering.'''
    
    k = 0
    
    while k < M.ncols():
        if M[i, k] < M[j, k]:
            return true
        elif M[i, k] > M[j, k]:
            return false
        k += 1
        
    return false

def sort_matrix(M):
    '''Sorts matrix rows lexicographically.'''
    
    for i in range(M.nrows()-1):
        for j in range(M.nrows()-i-1):
            if i_less_than_j(M, j, j+1):
                M = matrix(M.rows()[0:j] + M.rows()[j+1:j+2] + M.rows()[j:j+1] + M.rows()[j+2:])
                
    return M
            

class CAS(object):
    '''Defines class for Complete Access Structures.'''
    
    def _check_q_l(self, hw, dim, v):
        # Returns whether the access structure is Q_l (see
        # https://ia.cr/2018/467).

        # It does this by checking whether all vectors of
        # hamming-weight 'l' are contained in at least one maximally
        # unqualified set (i.e. in delta_plus)

        # hw = 'l' (as in Q_l)
        # dim = size of the vector left to populate in the recursive step (= number of parties)
        # v = Boolean vector encoding which parties are in the 'current'
        # set being checked (vector of dimension equal to number of
        # parties)

        if dim == -1:
            if hw == 0:
                i = 0
                found_all_parties = true
                while i < self.n and found_all_parties: # i.e. all parties so far are found
                    found_all_parties = false
                    for j in range(self.delta_plus.nrows()): # is_not_ql is true if at least one set contains i
                        if v[j] == 1:
                            found_all_parties = found_all_parties or self.delta_plus[j, i]
                    i += 1
                return found_all_parties
            else:
                return false

        v[dim] = 0
        if self._check_q_l(hw, dim - 1, v):
            return true

        v[dim] = 1
        if self._check_q_l(hw - 1, dim - 1, v):
            return true

        return false

    def _check_all_sets(self, dim, v):
        # Determines if each subset of parties is
        # qualified/unqualified and inserts in gamma/delta as
        # appropriate.

        # This ensures completeness of the access structure (i.e. that
        # every set is either qualified or unqualified

        if dim == -1:
            if self.delta[0, 0] != -1:
                j = 0
                found_subset = false
                while j < self.delta_plus.nrows() and not found_subset:
                    found_subset = true
                    i = 0
                    while i < self.n and found_subset:
                        if self.delta[j, i] < v[i]:
                            found_subset = false
                        i += 1
                    if found_subset and matrix(v) != self.delta[j, :]:
                        self.delta = matrix(self.delta.rows() + matrix(v).rows())
                    j += 1

            if self.gamma[0, 0] != -1:
                j = 0
                found_superset = false
                while j < self.gamma_minus.nrows() and not found_superset:
                    found_superset = true
                    i = 0
                    while i < self.n and found_superset:
                        if self.gamma[j, i] > v[i]:
                            found_superset = false
                        i += 1
                    if found_superset and matrix(v) != self.gamma[j, :]:
                        self.gamma = matrix(self.gamma.rows() + matrix(v).rows())
                    j += 1

            if not found_subset and not found_superset:
                if self._spec_as_unqualified and self.gamma[0, 0] != -1:
                    self.gamma = matrix(self.gamma.rows() + matrix(v).rows())
                elif self.delta[0, 0] != -1:
                    self.delta = matrix(self.delta.rows() + matrix(v).rows())
            return

        v[dim] = 0
        self._check_all_sets(dim - 1, v)

        v[dim] = 1
        self._check_all_sets(dim - 1, v)
    
        return 

    def _set_gamma_delta(self):
        # Calls _check_all_sets to populates gamma and delta.

        self._check_all_sets(self.n - 1, [0] * self.n)

        return

    def _set_delta_plus(self):
        # Puts maximal elements of self.delta in self.delta_plus.
        
        self.delta_plus = self.delta[0, :]

        for j in range(self.delta.nrows()):
            k = 0
            new_is_a_subset = false
            new_is_a_superset = false

            while k < self.delta_plus.nrows() and (not new_is_a_subset) and (not new_is_a_superset):
                i = 0
                new_is_a_superset = true
                while i < self.n and new_is_a_superset:
                    if self.delta[j, i] < self.delta_plus[k, i]:
                        new_is_a_superset = false
                    i += 1
                
                i = 0
                new_is_a_subset = true
                while i < self.n and new_is_a_subset:
                    if self.delta[j, i] > self.delta_plus[k, i]:
                        new_is_a_subset = false
                    i += 1
                k += 1
            k -= 1
                    
            if new_is_a_superset and (not new_is_a_subset):
                # Replace the first set of which new_superset is a
                # superset.  In the Python 'list' syntax below, we are
                # replacing the kth row with the new set.

                self.delta_plus = matrix(self.delta_plus.rows()[:k] + self.delta[j, :].rows() + self.delta_plus.rows()[(k+1):])

                # And remove subsets further down delta_plus
                k += 1
                while k < self.delta_plus.nrows():
                    i = 0
                    found_subset = true
                    while i < self.n and found_subset:
                        if self.delta[j, i] < self.delta_plus[k, i]:
                            found_subset = false
                        i += 1
                    if found_subset:
                        self.delta_plus = matrix(self.delta_plus.rows()[:k] + self.delta_plus.rows()[(k+1):])
                        k -= 1
                    k += 1

            if (not new_is_a_superset) and (not new_is_a_subset):
                self.delta_plus = matrix(self.delta_plus.rows() + self.delta[j, :].rows())
        return

    def _set_gamma_minus(self):
        # Puts minimal elements of self.gamma in self.gamma_plus.
        
        self.gamma_minus = self.gamma[0, :]

        for j in range(self.gamma.nrows()):
            k = 0
            new_is_a_subset = false
            new_is_a_superset = false
            while k < self.gamma_minus.nrows() and (not new_is_a_subset) and (not new_is_a_superset):
            #  Check if the current vector represents a superset of a currently alleged MUS
                i = 0
                new_is_a_superset = true
                while i < self.n and new_is_a_superset:
                    if self.gamma[j, i] < self.gamma_minus[k, i]:
                        new_is_a_superset = false
                    i += 1
                
                i = 0
                new_is_a_subset = true
                while i < self.n and new_is_a_subset:
                    if self.gamma[j, i] > self.gamma_minus[k, i]:
                        new_is_a_subset = false
                    i += 1

                k += 1
            k -= 1

            if (new_is_a_subset) and not (new_is_a_superset):
                # Replace the first set of which new_subset is a
                # subset.  In the Python 'list' syntax below, we are
                # replacing the kth row with the new set.
 
                self.gamma_minus = matrix(self.gamma_minus.rows()[:k] + self.gamma[j, :].rows() + self.gamma_minus.rows()[(k+1):])

                # And remove supersets further down gamma_minus
                k += 1
                while k < self.gamma_minus.nrows():
                    i = 0
                    found_superset = true
                    while i < self.n and found_superset:
                        if self.gamma[j, i] > self.gamma_minus[k, i]:
                            found_superset = false
                        i += 1
                    if found_superset:
                        self.gamma_minus = matrix(self.gamma_minus.rows()[:k] + self.gamma_minus.rows()[(k+1):])
                        k -= 1
                    k += 1

            if (not new_is_a_superset) and (not new_is_a_subset):
                self.gamma_minus = matrix(self.gamma_minus.rows() + self.gamma[j, :].rows())

        return

    def is_q(self, l):
        '''Determines if the access structure is Q_l (see https://ia.cr/2018/467).'''
        
        if l < self.n:
            return not self._check_q_l(l, self.delta_plus.nrows() - 1, list(0 for i in range(self.delta_plus.nrows())))
        
        return false
        
    def print_all(self):
        '''Prints all information about the access structure.'''
        
        print('Access structure info')
        print('=====================')
        display_matrix = matrix(QQ, matrix(QQ, list(i+1 for i in range(self.n))).rows() + self.delta.rows())
        print('Specified and deduced unqualified sets:\n', display_matrix)
        print
        display_matrix = matrix(QQ, matrix(QQ, list(i+1 for i in range(self.n))).rows() + self.delta_plus.rows())
        print('Derived maximally unqualified sets:\n', display_matrix)
        print
        display_matrix = matrix(QQ, matrix(QQ, list(i+1 for i in range(self.n))).rows() + self.gamma.rows())
        print('Deduced qualified sets:\n', display_matrix)
        print
        display_matrix = matrix(QQ, matrix(QQ, list(i+1 for i in range(self.n))).rows() + self.gamma_minus.rows())
        print('Derived minimally qualified sets:\n', display_matrix)
        print
        print('Number of parties:', self.n)
        print
        print('Number of maximally unqualified sets:', self.delta_plus.nrows())
        print
        print('Number of minimally qualified sets:', self.gamma_minus.nrows())
        print
        print('Is Q2?', self.is_q(2))

    def __init__(self, sets_to_use, unq):
        self.n = sets_to_use.ncols()
        self._spec_as_unqualified = unq

        if sets_to_use == matrix([-1] * self.n):
            if unq == 1:
                self.delta = matrix([-1] * self.n)
                self.delta_plus = matrix([-1] * self.n)
                self.gamma = matrix([0] * self.n)
                self.gamma_minus = matrix([0] * self.n)
                self._check_all_sets(self.n - 1, [0] * self.n)
            else:
                self.gamma = matrix([-1] * self.n)
                self.gamma_minus = matrix([-1] * self.n)
                self.delta = matrix([1] * self.n)
                self.delta_plus = matrix(QQ, [1] * self.n)
                self._check_all_sets(self.n - 1, [0] * self.n)
        elif unq == 1:
            self.delta = sets_to_use
            self.delta_plus = sets_to_use
            self.gamma = matrix([1] * self.n)
            self.gamma_minus = matrix([1] * self.n)
            self._check_all_sets(self.n - 1, [0] * self.n)
            if self.delta.nrows() == 2 ** self.n:
                self.gamma = matrix([-1] * self.n)
                self.gamma_minus = matrix([-1] * self.n)
            else:
                self.gamma = matrix(self.gamma.rows()[:self.gamma.nrows()])
        else:
            self.delta = matrix([0] * self.n)
            self.delta_plus = matrix([0] * self.n)
            self.gamma = sets_to_use
            self.gamma_minus = sets_to_use
            self._check_all_sets(self.n - 1, [0] * self.n)
            # Note that the next bit differs from the other case
            # because while the complete set of parties may still be
            # unqualified (although this is a pointless access
            # structure), the empty set of parties is always
            # unqualified.
            self.delta = matrix(self.delta.rows()[:self.delta.nrows()])

        self._set_delta_plus()
        self._set_gamma_minus()

        self.delta = sort_matrix(self.delta)
        self.delta_plus = sort_matrix(self.delta_plus)
        self.gamma = sort_matrix(self.gamma)
        self.gamma_minus = sort_matrix(self.gamma_minus)

class Threshold_CAS(CAS):
    '''Defines class for Threshold Complete Access Structures.'''

    def _generate_sets(self, hw, dim, v):
        if dim == -1:
            if hw == 0:
                self.delta = matrix(ZZ, matrix(ZZ, v).rows() + self.delta.rows())
            return

        v[dim] = 0
        self._generate_sets(hw, dim - 1, v)
        
        v[dim] = 1
        self._generate_sets(hw - 1, dim - 1, v)
        
        return

    def __init__(self, n, t):
        self.t = t
        self.n = n
        self.delta = matrix([0] * self.n)
        self._generate_sets(t, n - 1, [0] * n)
        self.delta = self.delta.delete_rows([self.delta.nrows()-1])
        CAS.__init__(self, self.delta, 1)
