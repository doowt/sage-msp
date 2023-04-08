load("cas.sage")

def hamming_weight(v):
    '''Returns the Hamming weight of a column vector.'''
    
    return int(v.norm(1))

class MSP(object):
    '''Defines class of Monotone Span Programs.'''

    def compute_all_msp_properties(self):
        # Computes and sets all properties of the MSP after
        # initialisation.  Takes a long time for large access
        # structures
        
        if self.target.ncols() != self.M.ncols():
            raise ValueError('Invalid target vector. Please set.')
        self._set_recomb_vectors()
        self._set_recomb_channels()
        self._set_error_matrix()
        self._set_list_of_reconst_matrices()
        return

    def _convert_msp(self, new_target):
        # Converts to an MSP that computes the same access structure
        # but has new_target as the target vector.

        if new_target.ncols() != self.M.ncols():
            text_to_output = 'The new target vector has dimension ' + str(new_target.ncols())
            test_to_output += ' whereas M has ' + str(self.M.ncols()) + ' columns.'
            raise ValueError(text_to_output)

        # First, find any part of the target which is non-zero
        k = 0
        found = false
        while k < self.target.ncols() and not found:
            if self.target[0, k] == 0:
                    k += 1
            else:
                found = true
        if not found:
            raise ValueError('The target vector is invalid: it must be a non-zero vector.')

        # Then perform row operations to modify the MSP matrix
        for i in range(new_target.ncols()):
            if new_target[0, i] == 0 and self.target[0, i] != 0:
                self.M[:, i] -= self.target[0, i] / self.target[0, k] * self.M[:, k]
            if new_target[0, i] != 0 and self.target[0, i] == 0:
                self.M[:, i] += new_target[0, i] / self.target[0, k] * self.M[:, k]
            if new_target[0, i] != 0 and self.target[0, i] != 0:
                self.M[:, i] = new_target[0, i] / self.target[0, i] * self.M[:, i]
            self.target[0, i] = new_target[0, i]
                
        return
    
    def _check_set_is_qualified(self, dim, v):
        # Used to determine an access structure from an MSP
        
        # Given a set of parties encoded as a vector v, check whether
        # or not it is qualified based on whether the target vector is
        # in the span of the rows owned by parties encoded in v.
        # Parameters are: dim = number of parties (decreases
        # recursively, starts at n) v = Boolean vector encoding which
        # parties are in the 'current' set being checked (vector of
        # dimension equal to number of parties)

        if dim == -1:
            matrix_to_test = self._sub_by_parties(v, self.M)
            if matrix(self.F, matrix_to_test.rows() + self.target.rows()).rank() == matrix_to_test.rank():
                self.AS.gamma = matrix(self.F, self.AS.gamma.rows() + matrix(self.F, v).rows())
        else:
            v[dim] = 0
            self._check_set_is_qualified(dim - 1, v)
            v[dim] = 1
            self._check_set_is_qualified(dim - 1, v)

    def compute_cas(self):
        '''Used to determine an access structure from an MSP.  Given an access structure, compute the complete access structure (where by assumption an undefined set is unqualified).'''

        self.AS.gamma = matrix([0] * self.AS.n)

        self._check_set_is_qualified(self.AS.n - 1, [0] * self.AS.n)

        if self.AS.gamma == matrix([0] * self.AS.n):
            self.AS.gamma = matrix([-1] * self.AS.n)
        else:
            self.AS.gamma = matrix(self.F, self.AS.gamma.rows()[1:])

        self.AS = CAS(self.AS.gamma, 0)
        
#        self.compute_all_msp_properties()

    def check_share_consistency(self, current_share):

        '''Checks that an input share vector is consistent by checking that for all shares corresponding to equal rows in M, the shares are equal.'''
        
        for i in range(self.M.nrows()):
            for j in range(i+1, self.M.nrows()):
                if self.M[i, :] == self.M[j, :] and current_share[i, 0] != current_share[j, 0]:
                    return false
        return true

    def check_share_correctness(self, current_share):
        '''Checks that an input share vector lies in the image of M.'''
        
        if self.N * current_share != matrix(self.F, [0] * self.N.nrows()).transpose():
            return false
        return true

    def _set_list_of_reconst_matrices(self):
        # Creates a reconstruction matrix for each party, to allow a
        # party to reconstruct the secret from the shares they observe
        # when a secret is 'opened'.
        
        self.list_of_reconst_matrices = []

        for i in range(self.AS.n):
            order = [0] * self.M.nrows()
            k = 0
            for j in range(self.M.nrows()):
                if self.comms_indicator[i, j] == 1:
                    order[k] = j
                    k += 1

            for j in range(self.M.nrows()):
                if self.comms_indicator[i, j] == 0:
                    order[k] = j
                    k += 1

            new_M = matrix(self.F, 1, self.M.ncols())
            for j in range(self.M.nrows()):
                new_M = matrix(self.F, new_M.rows() + matrix(self.F, self.M.rows()[order[j]]).rows())


            new_M = matrix(self.F, new_M.rows()[1:]).transpose().echelon_form().transpose()

            new_new_M = matrix(self.F, self.M.nrows(), self.M.ncols())

            for j in range(self.M.nrows()):
                for k in range(self.M.ncols()):
                    new_new_M[order[j], k] = new_M[j, k]

            self.list_of_reconst_matrices.append(new_new_M)

    def _check_v_in_span_A(self, v, A):
        # Checks to see if an input row vector v is in the rowspace of
        # A, and outputs a linear combination of the rows that
        # achieves v (if possible).

        aug_mat = matrix(self.F, A.nrows(), A.ncols() + A.nrows())

        for i in (range(A.nrows())):
            for j in range(A.ncols()):
                aug_mat[i, j] = A[i, j]
            aug_mat[i, A.ncols() + i] = 1

        aug_mat = aug_mat.echelon_form()
        
        inv_mat = matrix(self.F, A.nrows(), A.nrows())

        for i in range(A.nrows()):
            for j in range(A.nrows()):
                inv_mat[i, j] = aug_mat[i, j + A.ncols()]

        ech_recomb_vec = matrix(self.F, 1, aug_mat.nrows())

        # For each row
        for i in range(aug_mat.nrows()):
            # Find the first non-zero entry
            j = 0
            found = false
            while j < aug_mat.ncols() - inv_mat.ncols() and not found:
                if aug_mat[i, j] == 0:
                    j += 1
                else:
                    found = true
            # If it has one, add right multiple of the corresponding
            # part of the inverse matrix
            if found:
                ech_recomb_vec[0, i] = v[0, j]

        return ech_recomb_vec * inv_mat

    def _sub_by_parties(self, v, A):
        # Returns a submatrix of A whose rows are owned by parties
        # indicated by v.
        
        result = matrix(self.F, 1, A.ncols())
        for i in range(A.nrows()):
            if v[self.rowmap[i, 0]] == 1:
                result = matrix(self.F, result.rows() + matrix(self.F, A.row(i)).rows())
        return matrix(self.F, result.rows()[1:])

    def _sub_by_vector(self, v, A):
        # Returns a submatrix of A whose rows are indexed by v.
        
        result = matrix(self.F, 1, A.ncols())
        for i in range(A.nrows()):
            if v[i] == 1:
                result = matrix(self.F, result.rows() + matrix(self.F, A.row(i)).rows())
        return matrix(self.F, result.rows()[1:])

    def _sort_M_by_parties(self):
        # Sort M to make rowmap monotonically increasing.

        new_rowmap = matrix(self.F, self.M.nrows(), 1)
        new_M = matrix(self.F, self.M.nrows(), self.M.ncols())

        k = 0
        for i in range(self.AS.n):
            for j in range(self.M.nrows()):
                if self.rowmap[j, 0] == i:
                    new_rowmap[k, 0] = i
                    for l in range(self.M.ncols()):
                        new_M[k, l] = self.M[j, l]
                    k += 1

        self.M = new_M
        self.rowmap = new_rowmap

        self.compute_all_msp_properties()

    def _set_recomb_vectors(self):
        # Determines which shares each party should receive from the
        # other parties (see algorithm in appendix of
        # https://ia.cr/2018/467).  It returns a matrix in GF(2) of
        # dimension self.AS.n by self.M.nrows(), where (i,j) == 1 iff
        # party P_rowmap(j) sends P_i the j'th share (corresponding to
        # the j'th row in M).

        complete = matrix(self.F, 1, self.M.ncols())
        self.comms_indicator = matrix(self.F, self.AS.n, self.M.nrows())
        # First we will choose to look at the rows contributed by
        # minimally qualified sets of parties, because this will
        # usually give the desired result.  If this fails, we just
        # insert new rows from remaining parties until we have what we
        # need.

        # For load-balancing, we will search through gamma_minus such
        # that all the previously assigned sets have least preference
        # for subsequent parties.  This ensures, for example, that in
        # a (5,2)-threshold situation, party 1 receives from 2 and 3,
        # but party 2 receives from 3 and 4 instead of 1 and 3.
        set_used = [0] * self.AS.gamma_minus.nrows()
        for i in range(self.AS.n):
            # First we decide on the order to check through gamma_minus
            k = 0
            gamma_minus_order = [0] * self.AS.gamma_minus.nrows()
            for j in range(self.AS.gamma_minus.nrows()):
                if set_used[j] == 0:
                    gamma_minus_order[k] = j
                    k += 1
            for j in range(self.AS.gamma_minus.nrows()):
                if set_used[j] == 1:
                    gamma_minus_order[k] = j
                    k += 1

            # Then we use the order to look through gamma_minus for
            # the first set which contains the current party i
            mqs_index = 0
            found = false
            while mqs_index < self.AS.gamma_minus.nrows() and not found:
                if self.AS.gamma_minus[gamma_minus_order[mqs_index], i] == 1:
                    found = true
                    set_used[gamma_minus_order[mqs_index]] = 1
                    mqs_index -= 1
                mqs_index += 1

            # Now we will choose the order in which to insert new rows
            # to the new matrix (the one which is supposed to have
            # full rank).  After looking at all the rows owned by
            # parties in the minimally qualified set, we just insert
            # rows from any of the remaining parties until we have a
            # full-rank matrix.
            k = 1
            order_to_search = [0] * self.AS.n
            order_to_search[0] = i
            for j in range(self.AS.gamma_minus.ncols()):
                if self.AS.gamma_minus[gamma_minus_order[mqs_index], j] == 1 and j != i:
                    order_to_search[k] = j
                    k += 1

            for j in range(self.AS.gamma_minus.ncols()):
                if self.AS.gamma_minus[gamma_minus_order[mqs_index], j] == 0:
                    order_to_search[k] = j
                    k += 1

            # Now use order_to_search to find a matrix with full rank.
            current_matrix = matrix(self.F, 1, self.M.ncols())

            for j in range(self.AS.n): # i.e. length of order_to_search
                # Find the first row in M owned by the current party
                if current_matrix.rank() < self.M.ncols():
                    k = 0
                    found = false
                    while k < self.M.nrows() and not found:
                        if self.rowmap[k, 0] == order_to_search[j]:
                            found = true
                            k -= 1
                        k += 1

                    # Since the self.rowmap is 'block'ed into parties, just can
                    # increment by one each time instead of finding the
                    # next share owned by this party
                    while k < self.M.nrows() and found:
                        if self.rowmap[k, 0] != order_to_search[j]:
                            found = false
                        else:
                            if (matrix(self.F, current_matrix.rows() + matrix(self.F, self.M.row(k)).rows()).rank()) > current_matrix.rank():
                                current_matrix = matrix(self.F, current_matrix.rows() + matrix(self.F, self.M.row(k)).rows())
                                self.comms_indicator[i, k] = 1
                        k += 1

            current_recomb_vector = self._check_v_in_span_A(self.target, self._sub_by_vector(vector(self.comms_indicator.row(i)), self.M))

            complete = matrix(self.F, complete.rows() + current_recomb_vector.rows())

        self.recomb_vectors = matrix(self.F, complete.rows()[1:])
        self._set_recomb_channels()
        self._set_list_of_reconst_matrices()

    def _set_recomb_channels(self):
        # Sets the unidirectional channels required between parties to
        # open secrets.
        
        self.recomb_channels = matrix(self.F, self.AS.n, self.AS.n)
        
        for i in range(self.AS.n):
            for j in range(self.comms_indicator.ncols()):
                if self.comms_indicator[i, j] == 1 and self.rowmap[j, 0] != i:
                    self.recomb_channels[self.rowmap[j, 0], i] = 1

    def _set_error_matrix(self):
        # Sets self.N = self.M.left_kernel().

        aug_mat = matrix(self.F, self.M.nrows(), self.M.ncols() + self.M.nrows())

        for i in range(self.M.nrows()):
            for j in range(self.M.ncols()):
                aug_mat[i, j] = self.M[i, j]
            aug_mat[i, self.M.ncols() + i] = 1

        aug_mat = aug_mat.echelon_form()
        mat_rank = self.M.rank()

        self.N = matrix(self.F, aug_mat.nrows() - mat_rank, self.M.nrows())

        for i in range(aug_mat.nrows() - mat_rank):
            for j in range(self.M.nrows()):
                self.N[i, j] = aug_mat[mat_rank + i, self.M.ncols() + j]
                
        return 

    def _n_uni_directional_channels(self, mat):
        # Returns the number of unidirectional channels implied by a
        # matrix of channels.
        
        n_uni_channels = 0
        for j in range(mat.ncols()):
            n_uni_channels += int(hamming_weight(mat[:, j]))
        return n_uni_channels

    def _n_bi_directional_channels(self, mat):
        # Returns the number of bidirectional channels implied by the
        # input matrix of channels.
        
        n_bi_channels = 0
        for i in range(self.AS.n):
            for j in range(i + 1, self.AS.n):
                if mat[i, j] or mat[j, i]:
                    n_bi_channels += 1
        return n_bi_channels

    def print_all(self):
        '''Prints all information about the MSP.'''
        
        # self.AS.print_all()
        print('MSP info')
        print('========')
        print('Field:', self.F)
        print('(rowmap|MSP matrix):\n', self.rowmap.augment(self.M))
        print('Target:\n', self.target)
        print('Error-check matrix:\n', self.N)
        print('Recombination vectors:\n', self.recomb_vectors)
        print('Recombination channels:\n', self.recomb_channels)
        print('Total number of uni-directional channels used online:', self._n_uni_directional_channels(self.recomb_channels))
        print('Total number of bi-directional channels used online:', self._n_bi_directional_channels(self.recomb_channels))
       
    def __init__(self, AS, F):
        '''Initialises the MSP from an access structure AS and a field F.'''

        self.AS = AS
        self.F = F
        self.M = matrix(self.F)
        self.N = matrix(self.F)
        self.rowmap = matrix(self.F)
        self.target = matrix(self.F)
        self.recomb_vectors = matrix(self.F)
        self.recomb_channels = matrix(self.F)
        self.n_shares_held_by = [0] * self.AS.n
        
class Mult_MSP(MSP):
    '''Defines class of Multiplicative MSPs.'''
    
    def _compute_mult_as(self, tensored_matrix, tensored_target):
        # Computes and sets the access structure for the product of
        # two secrets.

        mult_msp = MSP(CAS(matrix(ZZ, [1] * self.AS.n), 1), self.F)

        mult_msp.M = matrix(self.F, tensored_matrix.rows())
        mult_msp.target = matrix(self.F, tensored_target.rows())
        mult_msp.rowmap = matrix(self.F, mult_msp.M.nrows(), 1)

        k = 0
        for i in range(self.AS.n):
            no_owned_by_i = 0
            for j in range(self.rowmap.nrows()):
                if self.rowmap[j, 0] == i:
                    no_owned_by_i += 1
            
            for j in range(no_owned_by_i ** 2):
                mult_msp.rowmap[k, 0] = i
                k += 1

        mult_msp.compute_cas()

        self.mult_AS = deepcopy(mult_msp.AS)
        
        return

    def _find_offline_mult_recomb_vector(self, dim, v):
        # Finds recombination vectors that allow parties to compute an
        # additive sum amongst possibly all parties of the product of
        # two secrets.
        
        if dim == -1:
            # Using the current v (which is just some element of
            # F_2^dim), search for a recombination vector in the
            # current_sub_matrix (which is the generator matrix
            # restricted to some fixed minimally qualified set) which
            # v defines (as a characteristic function for the rows)

            matrix_to_check = matrix(self.F, 1, self._current_sub_matrix.ncols())

            for i in range(self._current_sub_matrix.nrows()):
                if v[i] == 1:
                    matrix_to_check = matrix(self.F, matrix_to_check.rows() + matrix(self.F, self._current_sub_matrix.row(i)).rows())

            matrix_to_check = matrix(self.F, matrix_to_check.rows()[1:])

            recomb_sub_vector = self._check_v_in_span_A(self.target, matrix_to_check)

            if recomb_sub_vector * matrix_to_check == self.target:

                # Now make recomb_vector of the right dimension
                self.offline_mult_recomb_vector = matrix(self.F, 1, self._current_sub_matrix.nrows())
                k = 0
                for i in range(self._current_sub_matrix.nrows()):
                    if v[i] == 1:
                        self.offline_mult_recomb_vector[0, i] = recomb_sub_vector[0, k]
                        k += 1
        else:
            v[dim] = 0
            self._find_offline_mult_recomb_vector(dim - 1, v) 
            v[dim] = 1
            self._find_offline_mult_recomb_vector(dim - 1, v)


    def _convert_to_MMSP(self):
        # Converts an MSP to a multiplicative MSP, requiring that the
        # access structure be Q_2.
        
        ncols_in_dual = self.M.ncols()
        if self.AS.gamma_minus.nrows() < self.M.ncols():
            ncols_in_dual = self.AS.gamma_minus.nrows()

        # The theory is outlined by Cramer et
        # al. https://eprint.iacr.org/2000/037; it involves finding a
        # recombination vector for each minimally qualified set.

        # Roughly speaking, for each MQS, we determine a recombination
        # vector, which we insert as a column vector in the generator
        # matrix M (which we augment with some zeros to get linear
        # independence) where rowmap carries over in the natural way.

        # This method only works if the target vector is the all ones vector.
        self._convert_msp(matrix(self.F, [1] * self.M.ncols()))

        # Turn rowmap into two copies of itself
        new_rowmap = matrix(self.F, self.M.nrows() * 2, 1)
        for i in range(self.M.nrows()):
            new_rowmap[i, 0] = self.rowmap[i, 0]
            new_rowmap[i + self.M.nrows(), 0] = self.rowmap[i, 0]

        # M must be widened by the number of MQSs (i.e. AS.gamma_minus.nrows())
        # and doubled in height.  We create a new matrix new_M and
        # copy in the data.

        new_M = matrix(self.F, self.M.nrows()*2, self.M.ncols() + ncols_in_dual - 1)
        for i in range(self.M.nrows()):
            for j in range(self.M.ncols()):
                new_M[i, j] = self.M[i, j]

        # Now extract the first d qualified sets from the access structure
        for i in range(ncols_in_dual):

            # And get the submatrix of M indexed by the parties in this set...
            self._current_sub_matrix = self._sub_by_parties(list(self.AS.gamma_minus.row(i)), self.M)

            # ...in order to find a (projected) recombination vector
            # within this submatrix
            self.offline_mult_recomb_vector = matrix(self.F, 1, self._current_sub_matrix.nrows())

            self._find_offline_mult_recomb_vector(self._current_sub_matrix.nrows() - 1, [0] * self._current_sub_matrix.nrows())

            # Lift the vector into the right dimension and insert into
            # M in the right position

            k = 0
            for j in range(self.M.nrows(), self.M.nrows()*2):
                if self.AS.gamma_minus[i, new_rowmap[j, 0]] == 1:
                    if i == 0:
                        new_M[j, 0] = self.offline_mult_recomb_vector[0, k]
                    else:
                        new_M[j, self.M.ncols() + i - 1] = self.offline_mult_recomb_vector[0, k]
                    k += 1

        # We have to do some transformations of M to get the correct
        # output. Since this is just column operations on M, M*, and
        # the joint MSP, this does not change the access structure
        # they compute.

        for i in range(self.M.nrows()):
            for j in range(self.M.ncols(), new_M.ncols()):
                new_M[i, j] += new_M[i, 0]

        for i in range(self.M.nrows(), new_M.nrows()):
            for j in range(1, self.M.ncols()):
                new_M[i, j] += new_M[i, 0]

        self.target = matrix(self.F, [1] * new_M.ncols())

        self.M = new_M
        self.rowmap = new_rowmap
        self._sort_M_by_parties()
        self.compute_all_msp_properties()

        return
    
    def _is_multiplicative(self, current_MSP):
        # Returns whether or not an MSP is multiplicative.
        
        complete = matrix(self.F, 1, current_MSP.M.ncols()*current_MSP.M.ncols())
        for i in range(current_MSP.AS.n):

            current_party = [0] * current_MSP.AS.n
            current_party[i] = 1
            current_matrix = self._sub_by_parties(current_party, self.M)
            new_row = matrix(self.F, 1, current_matrix.ncols() * current_matrix.ncols())
            tensor_matrix = matrix(self.F, 1, current_matrix.ncols() * current_matrix.ncols())
            
            # Compute tensors the hard way
            for j in range(current_matrix.nrows()):
                for k in range(current_matrix.nrows()):
                    for l in range(current_matrix.ncols()):
                        for m in range(current_matrix.ncols()):
                            new_row[0, l*current_matrix.ncols() + m] = current_matrix[j, l] * current_matrix[k, m]
                    tensor_matrix = matrix(self.F, tensor_matrix.rows() + new_row.rows())
            complete = matrix(self.F, complete.rows() + tensor_matrix.rows()[1:])

        # Remove the first row since it is an artefact of
        # initialisation and is the zero row
        complete = matrix(self.F, complete.rows()[1:])
        ttarget = matrix(self.F, 1, current_MSP.target.ncols() * current_MSP.target.ncols())

        for i in range(current_MSP.target.ncols()):
            for j in range(current_MSP.target.ncols()):
                ttarget[0, i*current_MSP.target.ncols() + j] = current_MSP.target[0, i] * current_MSP.target[0, j]
                # TODO: This might be the wrong way around -- check!
                # Note that it doesn't make a difference for the two
                # most common target vectors

        self.offline_mult_recomb_vector = self._check_v_in_span_A(ttarget, complete)

        self.list_of_mult_matrices = []

        posn_in_offline_mult_recomb_vector = 0
        for i in range(self.AS.n):
            
            # First count the number of shares this party holds
            no_of_shares = 0
            for j in range(self.rowmap.nrows()):
                if self.rowmap[j, 0] == i:
                    no_of_shares += 1

            # Since we ran self._sort_M_by_parties(), the first
            # no_of_shares^2 entries of the vector
            # self.offline_mult_recomb_vector are the coefficients
            # corresponding to local products this party computes.
            current_matrix = matrix(self.F, no_of_shares, no_of_shares)
            current_index = 0
            for j in range(no_of_shares):
                for k in range(no_of_shares):
                    current_matrix[j, k] = self.offline_mult_recomb_vector[0, posn_in_offline_mult_recomb_vector + current_index]
                    current_index += 1
            self.list_of_mult_matrices.append(current_matrix)
            posn_in_offline_mult_recomb_vector += no_of_shares*no_of_shares

        if self.offline_mult_recomb_vector * complete == ttarget:
            self._compute_mult_as(complete, ttarget)
        else:
            return false

        return true

    def __init__(self, current_MSP):
        '''Returns the input MSP if already multiplicative, and if not then computes and returns a multiplicative MSP.'''

        # Initialise the output as whatever was given as input
        if not current_MSP.AS.is_q(2):
            raise ValueError('Access structure not Q2 -- cannot be made multiplicative.')

        MSP.__init__(self,current_MSP.AS, current_MSP.F)

        self.M = current_MSP.M
        self.N = current_MSP.N
        self.rowmap = current_MSP.rowmap
        self.target = current_MSP.target
        

        self._sort_M_by_parties()

        if not self._is_multiplicative(current_MSP):
            self._convert_to_MMSP()
            self._is_multiplicative(self) # To recompute self.offline_mult_recomb_vector

    def print_all(self):
        '''Prints all information about the Multiplicative MSP.'''
        
        MSP.print_all(self)
        print()
        print('m-MSP-specific info')
        print('========================')
        print('Matrices for computing the additive sharing of product:')
        for i in range(self.AS.n):
            print(i, ": \n", self.list_of_mult_matrices[i])

class Replicated_MSP(MSP):
    '''Defines a class for replicated secret-sharing (Ito, M., Saito, A. and Nishizeki, T. (1987) Secret Sharing Scheme Realizing General Access Structure. Proceedings of IEEE Globecom 87, pp99-102.).'''
    
    def _set_already_assigned(self, indexes, i):
        # Returns (party i has already been assigned a set).
        
        for j in range(i):
            if indexes[i] == indexes[j] and i != j:
                return true
        return false

    def _number_sets_assigned(self, indexes, limit):
        # Returns number of sets so far assigned to parties (all must
        # be assigned at least once for opening secrets).
        
        sets_assigned = [0] * self.AS.delta_plus.nrows()
        for i in range(self.AS.n): # limit+1):
            if indexes[i, 0] > -1 and indexes[i, 0] < self.AS.delta_plus.nrows():
                sets_assigned[indexes[i, 0]] = 1
        n_assigned = 0
        for i in range(self.AS.delta_plus.nrows()):
            if sets_assigned[i] == 1:
                n_assigned += 1
        return n_assigned

    def _find_one_assignment(self, indices, i):
        # Determines if an assignment of parties to sets (shares) is
        # possible, and if so, sets self.parties_to_sets to be this
        # map
        
        # -2 means party not assigned a set (yet)
        # -1 means null

        if self._number_sets_assigned(indices, self.AS.n-1) == self.AS.delta_plus.nrows():
            self.parties_to_sets = copy(indices)
            return true
        elif i == self.AS.n:
            if self._number_sets_assigned(indices, self.AS.n-1) >= self._number_sets_assigned(self.parties_to_sets, self.AS.n-1):
                self.parties_to_sets = copy(indices)
        else:
            if i < self.AS.n - 1:
                next_non_null = i + 1
                while indices[next_non_null, 0] == -1 and next_non_null < self.AS.n:
                    next_non_null += 1
            else:
                next_non_null = self.AS.n

            if i > 0:
                previous_non_null = i - 1
                while indices[previous_non_null, 0] == -1 and previous_non_null > -1:
                    previous_non_null -= 1
            else:
                previous_non_null = -1

            found_set = false

            if indices[i, 0] == -2:
                indices[i, 0] = 0

            while indices[i, 0] < self.AS.delta_plus.nrows() and not found_set:
                if ((self.AS.delta_plus[indices[i, 0], i] == 0) and (not self._set_already_assigned(indices, i))):
                    found_set = true
                else:
                    indices[i, 0] += 1

            if found_set:
                self._find_one_assignment(indices, next_non_null)
            elif i > 0:
                temp_indices = copy(indices)
                temp_indices[i, 0] = -1
                if not self._find_one_assignment(temp_indices, next_non_null):
                    indices[i, 0] = -2
                    temp_indices = copy(indices)
                    if previous_non_null > -1:
                        temp_indices[previous_non_null, 0] += 1
                        self._find_one_assignment(temp_indices, previous_non_null)
                    else:
                        temp_indices[i, 0] = -1
                        self._find_one_assignment(temp_indices, next_non_null)
            else:
                temp_indices = copy(indices)
                temp_indices[i, 0] = -1
                self._find_one_assignment(temp_indices, next_non_null)
                
        return false
            
    def _assign_all_remaining(self, assignment):
        # Assign any unassigned sets to arbitrary parties so that each
        # share will be opened by at least one party

        reverse_assignment = matrix(self.F, [-2] * self.AS.delta_plus.nrows()).transpose()
        for i in range(self.AS.n):
            if assignment[i, 0] > -1 and assignment[i, 0] < self.AS.delta_plus.nrows():
                reverse_assignment[assignment[i, 0], 0] = i
                
        for j in range(self.AS.delta_plus.nrows()):
            if reverse_assignment[j, 0] == -2:
                i = 0
                while self.AS.delta_plus[j, i] != 0 and i < self.AS.n - 1:
                    i += 1
                if i == self.AS.n - 1 and self.AS.delta_plus[j, i] != 0:
                    reverse_assignment[j, 0] = -1
                else:
                    reverse_assignment[j, 0] = i

        return reverse_assignment

    def _set_sets_to_parties(self):
        # Sets the map self.sets_to_parties to encode which party must
        # open which share.
        
        # Find one map of parties to sets, and set
        # self.parties_to_sets if it exists
        self.parties_to_sets = matrix(self.F, [-2] * self.AS.n).transpose()
        self._find_one_assignment(matrix(self.F, [-2] * self.AS.n).transpose(), 0)

        # Now reverse the map to find self.sets_to_parties
        self.sets_to_parties = matrix(self.F, [-2] * self.AS.delta_plus.nrows()).transpose()
        for i in range(self.AS.n):
            for j in range(self.AS.n):
                if self.parties_to_sets[i] == j:
                    self.sets_to_parties[j] = i

        self.sets_to_parties = self._assign_all_remaining(self.parties_to_sets)
        return

    def print_all(self):
        '''Prints all information about the Replicated MSP.'''
        
        MSP.print_all(self)
        print()
        print('Replicated-specific info')
        print('========================')
        print('Map of sets to parties (to define which share is opened by which party):\n', self.sets_to_parties)
        print('Map of parties to sets:\n', self.parties_to_sets)

    def __init__(self, AS, F):
        '''Initialises instance of replicated secret-sharing from an access structure AS and a field F.'''

        MSP.__init__(self, AS, F)
        
        shares = matrix(self.F, self.AS.delta_plus.nrows(), self.AS.delta_plus.nrows())
        for i in range(self.AS.delta_plus.nrows()):
            # shares[0, i] = -1
            shares[i, i] = 1

        # Set the target vector
        self.target = matrix(self.F, 1, self.AS.delta_plus.nrows())
        for i in range(self.AS.delta_plus.nrows()):
            self.target[0, i] = 1

        # Set M and the rowmap
        m = 0
        for i in range(self.AS.delta_plus.nrows()):
            m += self.AS.n - int(hamming_weight(self.AS.delta_plus[i, :].transpose()))

        self.M = matrix(self.F, m, self.AS.delta_plus.nrows())
        self.rowmap = matrix(self.F, m,1)
        k = 0
        for i in range(self.AS.n):
            for j in range(self.AS.delta_plus.nrows()):
                if self.AS.delta_plus[j, i] == 0:
                    self.M[k, :] = shares[j, :]
                    self.rowmap[k, 0] = i
                    k += 1

        # Set the rest
        self._sort_M_by_parties()
        self._set_sets_to_parties()
        self.compute_all_msp_properties()

        self.n_shares_held_by = [0] * self.AS.n
        for current_party in range(self.AS.n):
            for j in range(self.rowmap.nrows()):
                if self.rowmap[j, 0] == current_party:
                    self.n_shares_held_by[current_party] += 1
                    

class Shamir_MSP(MSP):
    '''Defines a class for Shamir's secret-sharing (Shamir, A. How to Share a Secret, 1979).'''
    
    def __init__(self, AS, F):
        '''Initialises instance of Shamir's secret-sharing from an access structure AS and a field F.'''

        MSP.__init__(self, AS, F)

        # Set the target vector
        self.target = matrix(self.F, 1, hamming_weight(self.AS.delta_plus[0, :].transpose())+1)
        self.target[0, 0] = 1

        # Set M and the rowmap
        v = list(i+1 for i in range(self.AS.n))

        # TODO: Allow non-canonical definition of shares (choose
        # the field elements) or use a primitive root raised to
        # various powers (which implies efficient decoding)
        self.M = matrix(self.F, self.AS.n, hamming_weight(self.AS.delta_plus[0, :].transpose())+1)
        for i in range(self.AS.n):
            for j in range(hamming_weight(self.AS.delta_plus[0, :].transpose())+1):
                self.M[i, j] = (i+1) ** j

        self.rowmap = matrix(self.F, list(i for i in range(self.AS.n))).transpose()

        self.n_shares_held_by = [1] * self.AS.n

        # Set the rest
        self._sort_M_by_parties()
        self.compute_all_msp_properties()
