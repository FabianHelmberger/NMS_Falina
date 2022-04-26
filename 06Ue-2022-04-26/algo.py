def Gershgorin(A):
    dim = np.shape(A)[0]
    Radius = []

    out = []

    # loop over diagonal elemts
    for j in range(dim):
        tuple = []

        # loop oder cols and rows
        for dir in ['col', 'row']:
            if dir == 'row':
                A = A.T

            diff = 0
            for i in range(dim):

                # don't sum diagonals
                if i != j:
                    diff = diff + np.abs(A[j,i])

            tuple.append([diff])
        out.append(tuple)
    return out