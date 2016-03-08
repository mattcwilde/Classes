def sort(a, lo, hi, aux):
    # base case of recursion to end recursion
    if (hi <= lo):
        return

    #      low        mid        hi
    # []   []   []   []    []    []    []
    mid = lo + (hi - lo)/2

    # calling a function in a definition of that function
    # this is what recursion is

    # DIVIDE and CONQUER!!!!!
    # sort left half

    # if the parameters were equal to the initial parameters it would run forever
    sort(a, lo, mid, aux) 

    # sort right half
    sort(a, mid+1, hi, aux) 

    #now merge the two sorted halfs
    merge(a, lo, mid, hi, aux)


def merge(a, lo, mid, hi, aux):
    i = lo
    j = mid + 1

    # make a copy
    #aux = [0]*len(a)
    aux[lo:hi+1] = a[lo:hi+1]

    k = lo
    while (k <= hi):
        while((i <= mid) and (j <= hi)):

            if(aux[j] < aux[i]):
                #pick the smaller one
                a[k] = aux[j]
                k = k + 1
                j = j + 1
            else:
                a[k] = aux[i]
                i = i + 1 
                k = k + 1 

        # Now either i > mid or j > hi
        # so only one of the below loops will run
        while(i <= mid):
            a[k] = aux[i]
            i = i + 1
            k = k + 1 
        while(j <= hi):
            a[k] = aux[j]
            j = j + 1
            k = k + 1
    return a  

def mergesort(a):
    aux = [0]*len(a)
    # overloading funciton ok since takes a different amount of input
    # signature of function is different so computer names them differently
    sort(a, 0, len(a)-1, aux)
    return a