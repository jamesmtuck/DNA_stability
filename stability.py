# Contributed by James Tuck (jtuck@ncsu.edu) 

import sys
import csv
import math
from math import factorial

import matplotlib.pyplot as plt
import gmpy2 as g
from gmpy2 import mpfr

# Tune precision of gmpy2 module
g.get_context().precision = 1000

def dumpXY(name, XY, labels):
    """ Write a csv file called name+".csv" using
        data in XY organized as tuples (X0,Y0), (X1,Y1), 
        ... (Xi, Yi), where each Xi and Yi are vectors of a 
        common length, n.
        Each line of the file expands these sets:
            X0[0], Y0[1], X0[1], Y0[1],...,X0[j],Y0[j],...X0[n-1],Y0[n-1]
            X1[0], Y1[1], X1[1], Y1[1],...,X1[j],Y1[j],...X1[n-1],Y1[n-1]
            ...
            Xi[0], Yi[1], Xi[1], Yi[1],...,Xi[j],Yi[j],...Xi[n-1],Yi[n-1]
            ...
    """
    with open(name+".csv","w") as csvfile:
        wr = csv.writer(csvfile)
        row = []
        for l in labels:
            row.append("x")
            row.append(l)
        wr.writerow(row)
        row = []
        l = 0
        for x,y in XY:
            l = max(l, len(x))
        for i in range(l):
            row = []
            for x,y in XY:
                if i < len(x):
                    row.append(x[i])
                    row.append(y[i])
                else:
                    row.append(' ')
                    row.append(' ')
            wr.writerow(row)

    
def dump(name,X,X_label,Ys,Ys_labels):
    """ Write a csv file called name+".csv" using
        X as independent variable and Ys as a set of dependent variables.
        xlabel is the label for X, and ylabels for Y.
        CSV file contains lines like this, where m is length of X and Y and n is length
        of each member of Y:
        X[0],Y[0][0],Y[0][1],Y[0][2],...,Y[0][n-1]
        X[1],Y[1][0],Y[1][1],Y[1][2],...,Y[1][n-1]
        ...
        X[i],Y[i][0],Y[i][1],Y[i][2],...,Y[i][n-1]
        ...
        X[m-1],Y[m-1][0],Y[m-1][1],Y[m-1][2],...,Y[m-1][n-1]

    """
    with open(name+'.csv', 'w') as csvfile:
        wr = csv.writer(csvfile)
        wr.writerow([X_label]+Ys_labels)
        for i,x in enumerate(X):
            row = [x]
            for y in Ys:
                row.append("{:1.5e}".format(y[i]))
            wr.writerow(row) 


def plot(X,Ys,labels,xlabel="",ylabel="",title=""):
    """ Helper function to quickly plot Y vs X with corresponding labels. """
    for Y,label in zip(Ys,labels):
        plt.plot(X,Y,label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.suptitle(title)
    plt.legend()
    plt.show()

def g_pow(base,e):
    """ g_pow(base,e): gmpy2 helper function for exponentiation
        with proper special case handling for base==0.
    """
    #print (base,e)
    if base==mpfr("0.0"):
        if e==mpfr(0):
            return mpfr("1")
        else:
            return mpfr("0")
    return g.exp(g.mul(g.log(mpfr(base)),mpfr(e)))


def binomial_pmf(q, N, k):
    """ Calculate binomial probability mass function. """
    # Useful as sanity check that gmpy2 is providing sufficient precision.
    # g.bincoef is essential for ensuring precision.
    tmp = g.mul(g_pow(q,k),g.mul(g.bincoef(N,k),g_pow(1-q,N-k)))
    return tmp

def binomial_cdf(q, N, k):
    """ Calculate binomial cumulative distribution function. """
    # Useful as sanity check that gmpy2 is providing sufficient precision.
    # g.bincoef is essential for ensuring precision.
    tmp_list = [mpfr("0")]
    for i in range(0,k+1):
        tt1 = g.mul(g_pow(q,i),g.mul(g.bincoef(N,i),g_pow(1-q,N-i)))
        tmp_list.append( tt1 )   
    tmp1 = g.fsum(tmp_list)
    return tmp1

def plot_binomial_pmf(N):
    X = [ _ for _ in range(0,N+1,1) ]
    Y = [
        [ binomial_pmf(0.01,N,_) for _ in X ],
        [ binomial_pmf(0.25,N,_) for _ in X ],
        [ binomial_pmf(0.5,N,_) for _ in X ],
        [ binomial_pmf(0.7,N,_) for _ in X ],
        [ binomial_pmf(0.9,N,_) for _ in X ]
         ]
    #print (Y)
    plot(X,Y,[ "q={} N={}".format(qi,N) for qi in [0.01,0.25,0.5,0.7,0.9] ],
         "Random Variable", "Probability", "Binomial PMF with N={} and varying q".format(N))


def plot_binomial_cdf(N):
    X = [ _ for _ in range(0,N+1,1) ]
    Y = [
        [ binomial_cdf(0.01,N,_) for _ in X ],
        [ binomial_cdf(0.25,N,_) for _ in X ],
        [ binomial_cdf(0.5,N,_) for _ in X ],
        [ binomial_cdf(0.7,N,_) for _ in X ],
        [ binomial_cdf(0.9,N,_) for _ in X ]
         ]
    #print (Y)
    plot(X,Y,[ "q={} N={}".format(qi,N) for qi in [0.01,0.25,0.5,0.7,0.9] ], "Random Variable", "Probability",
         "Binomial CDF with N={} and varying q".format(N))


def plot_binomials(N):
    """ Plot binomial PMF and CDF for value of N """
    plot_binomial_cdf(N)
    plot_binomial_pmf(N)

class Codeword:
    """ This class models a codeword very abstractly. """
    def __init__(self,length,num,es=1e-3):
        """
           Its members are:
              num = number of codewords
              length = length of codeword in nucleotides, assumes a block code
              es = probability of single nucleotide in error
        """
        self.length = length
        self.num = num
        self.es = es

    def P_error(self,k=1):
        """ Compute the likelihood of a codeword being in error based on es and length """
        return g.sub(mpfr("1"),binomial_cdf(self.es,self.length,k-1))

class RSCode:
    """ RSCode models a single Reed-Solomon Code, either inner or outer. """
    
    """ memoize RSCode objects so that we avoid unnecessary repeated calculations """
    rs_code_objects = {}

    @staticmethod
    def create(n,k,d,es=1e-3,ee=1e-3):
        """  Build an object that represents an RS[n,k,d] code, and specify the probability
             of a symbol error (es) or erasure (ee).         
        """
        if (n,k,d,es,ee) in RSCode.rs_code_objects:
            return RSCode.rs_code_objects[(n,k,d,es,ee)]
        else:
            rs = RSCode(n,k,d,es,ee)
            RSCode.rs_code_objects[(n,k,d,es,ee)] = rs
            return rs

    def __init__(self,n,k,d,es=1e-3,ee=1e-3):
        """ Constructor for n [n,k,d] RSCode with symbol error rate (es) and erasure
            error rate (ee).
        """
        self.q = 4
        self.n = n
        self.k = k
        self.d = d        
        self.t = int((d-1)/2)
        self.symbol_err_rate = es
        self.erasure_err_rate = ee
        self.result = mpfr("0")
        self.has_result = False
        #print (n,k,d,es,ee)

    def label(self):
        return "RS[{},{},{}] psym={:1.2e} perasure={:1.2e}".format(self.n,self.k,self.d,self.symbol_err_rate, self.erasure_err_rate)
    
    def R(self):
        return self.k / self.n

    def R_index(self, M):
        return (self.k - math.log(M,256))/float(self.n)
    
    def P_symbol_error(self):
        return g.sub(mpfr("1"),binomial_cdf(self.symbol_err_rate,self.n,self.t))
        
    def P_random_errors(self,i):
        return binomial_pmf(self.symbol_err_rate,self.n,i)
        
    def P_random_erasures(self,x):
        #print ("random_erasures",x,g.sub(1,binomial_cdf(self.erasure_err_rate,self.n,x)))
        return g.sub(mpfr("1"),binomial_cdf(self.erasure_err_rate,self.n,x))
    
    def P_erasure_error(self):
        tmp = mpfr('0')
        for i in range(self.t+1):
            tmp = g.add(tmp, g.mul(self.P_random_errors(i),self.P_random_erasures(2*(self.t-i))))
        return tmp
    
    def P_result(self):
        if not self.has_result:
            self.result = g.add(self.P_symbol_error(),self.P_erasure_error())
            self.has_result = False
        return self.result

    #def P_erasure_res(self):
    #    return self.P_random_erasures(self.d-1)
        
    #def P_ue(self):
    #    return g.sub(mpfr("1"),binomial_cdf(self.symbol_err_rate,self.n,self.t))
    

class RSInnerOuter:
    """ RSInnerOuter contains an inner and outer model RSCode object and is able
        to analyze their combined properties, namely code Rate and decoding error
        probability.
    """
    def __init__(self,n_inner,k_inner,d_inner,n_outer,k_outer,d_outer,p_sub,p_strand_loss):
        self.rs_inner = RSCode.create(n_inner,k_inner,d_inner,p_sub,mpfr("1e-15"))
        #print ("[{},{},{}] Error probability inner={:1.2e}".format(n_inner,k_inner,d_inner,self.rs_inner.P_result()))
        self.rs_outer = RSCode.create(n_outer,k_outer,d_outer,self.rs_inner.P_result(),mpfr(p_strand_loss))
        #print ("[{},{},{}] Error probability outer={:1.2e}".format(n_outer,k_outer,d_outer,self.rs_outer.P_result()))
        
    def P_result(self):
        return self.rs_outer.P_result()

    def R_index(self, M):
        #print ("outer_k={} inner_k={} logM={}".format(self.rs_outer.k,self.rs_inner.k,math.log(M,256)))
        return (self.rs_outer.k)*(self.rs_inner.k - math.ceil(math.log(M, 256))) / (self.rs_inner.n*self.rs_outer.n)
    
    def R_raw(self):
        return (self.rs_inner.k * self.rs_outer.k)/(self.rs_inner.n*self.rs_inner.n)

    def getLabel(self):
        return "[{}*{},{}*{}]".format(self.rs_outer.n,self.rs_inner.n,self.rs_outer.k,self.rs_inner,k)

    # @staticmethod
    # def get(L,d_inner,d_outer,p_sym,p_strand):
    #     print ("get args:",L,d_inner,d_outer,p_sym,p_strand)
    #     cw = Codeword(4,256,p_sym)
    #     n_inner = int(L/4)
    #     k_inner = n_inner - d_inner
    #     n_outer = 255
    #     k_outer = n_outer - d_outer
    #     return RSInnerOuter(n_inner,k_inner,d_inner,n_outer,k_outer,d_outer,p_sym,p_strand)    
    
    # @staticmethod
    # def get50(p_sub,p_strand_loss):
    #     expected_strand_loss = int(255*p_strand_loss*(1-p_strand_loss)*2+1)
    #     cw = Codeword(4,256,p_sub)
    #     L = 50
    #     n_inner = int(L/4)
    #     d_inner = 2*int(n_inner*cw.P_error()*(1-cw.P_error))+1
    #     k_inner = n_inner - d_inner
    #     n_outer = 255
    #     d_outer = expected_strand_loss
    #     k_outer = n_outer - d_outer
    #     return RSInnerOuter(n_inner,k_inner,d_inner,n_outer,k_outer,d_outer,p_sub,p_strand_loss)

    # @staticmethod
    # def get200(p_sub,p_strand_loss):
    #     expected_strand_loss = int(255*p_strand_loss*(1-p_strand_loss)*2+1)
    #     cw = Codeword(4,256,p_sub)
    #     L = 200
    #     n_inner = int(L/4)
    #     d_inner = 2*int(n_inner*cw.P_error()*(1-cw.P_error))+1
    #     k_inner = n_inner - d_inner
    #     n_outer = 255
    #     d_outer = expected_strand_loss
    #     k_outer = n_outer - d_outer
    #     return RSInnerOuter(n_inner,k_inner,d_inner,n_outer,k_outer,d_outer,p_sub,p_strand_loss)

    # @staticmethod
    # def get1000(p_sub,p_strand_loss):
    #     expected_strand_loss = int(255*p_strand_loss*(1-p_strand_loss)*2+1)
    #     cw = Codeword(4,256,p_sub)
    #     L = 1000
    #     n_inner = int(L/4)
    #     d_inner = 2*int(n_inner*cw.P_error()*(1-cw.P_error))+1
    #     k_inner = n_inner - d_inner
    #     n_outer = 255
    #     d_outer = expected_strand_loss
    #     k_outer = n_outer - d_outer
    #     return RSInnerOuter(n_inner,k_inner,d_inner,n_outer,k_outer,d_outer,p_sub,p_strand_loss)

      
def result(L,cw_size=4,cw_er=1e-3,cw_es=1e-3,dropout=1e-3):
    n = int(L/cw_size)
    errs = []
    D = []
    for d in range(3,n,2):
        k = n - d
        D.append(d)
        rs_inner = RSCode(n*cw_size,k*cw_size,d,cw_er,cw_es)
        rs_outer = RSCode(255,255-31,31,rs_inner.P_result(),dropout)
        errs.append(g.log(rs_outer.P_result()))
    #print (zip(D,errs))
    plot (D,[errs],["p"])    
    return


def rs_sweep_length(n=255,k=[235],Length=[_ for _ in range(50,201,25)],p_se=[g_pow(10,_/10.0) for _ in range(-80,-10,1)],copies=[1],PL=lambda L, c: g_pow(g.mul(mpfr(1e-3),mpfr(L)),mpfr(c)), filename="sweep_length" ):

    #print (n,k,Length,p_se)    
    #X = [ math.log(x,10) for x in p_se ]
    X = Length
    #print (X)
    Y = []
    Label = []
    for ki in k:
        #nol = []
        for c in copies:
            for p in p_se:
                y = []
                Label.append( "RS(N={},k={}) perr={:1.1e} c={}".format(n,ki,float(p),c) )
                for l in Length:
                    #cw_size = 4
                    #cw_strand = int(l/cw_size)
                    #d = 3                
                    #rs_inner = RSCode(cw_strand,cw_strand-d,d,mpfr(p),mpfr(0.001))
                    rs_code = RSCode(n,ki,n-ki,mpfr(p),PL(mpfr(l),c))
                    y.append(g.log10(rs_code.P_result()))
                Y.append(y)
                #print (y)                        

    fig, ax = plt.subplots()
    #ax.set_xlim(-1,-8)
    #ax.set_ylim(-20,0)
    plt.ylabel("Probabilty Decoding Error - Log Scale (10^y)")
    plt.xlabel("Strand Length (nt)")
    #ax.grid(True)
    for y,label in zip(Y,Label):
        plt.plot(X,y,label=label)
    plt.legend()
    plt.show()
    dump(filename,X,"Length",Y,Label)


def get_inner_at_target(L,index,P_target,p_sym):
    """ Search for an inner RSCode that is at least as good as P_target
        assuming L and symbol error rate of p_sym. 
    """
    cw = Codeword(4,256,p_sym)
    n = int(L/4)
    inner = None
    guess = int(n*p_sym)
    if guess % 2 == 0:
        guess += 1
    for d in range(guess,n-3,2):
        if n-d-index < 1:
            continue
        inner = RSCode.create(n,n-d,d,es=mpfr(cw.P_error()),ee=mpfr(1e-15))
        if inner.P_result() <= P_target:
            return inner
    return inner

def get_outer_at_target(P_target,p_sym,p_loss):
    """ Search for an outer RSCode by varying distance that is at least 
        as good as P_target assuming L=255 and symbol error rate of p_sym 
        and strand loss rate of p_loss 
    """
    n=255
    outer = None
    guess = int(n*p_loss)
    if guess % 2 == 0:
        guess += 1
    for d in range(guess,n-3,2):
        outer = RSCode.create(n,n-d,d,p_sym,p_loss)
        if outer.P_result() < P_target:
            return outer
    return outer
    
    
def get_rs_at_res(M=1e9,L=200,P_res=mpfr(1e-11),p_sym=mpfr(1e-3),p_loss=mpfr(1e-3)):
    cw = Codeword(4,256,p_sym)
    index = math.log(M,256)
    inner = get_inner_at_target(L,index,P_res,p_sym)
    if inner==None:
        return None
    frac = g.factorial(int((inner.d-1)/2))
    outer = get_outer_at_target(P_res,g.div(inner.P_result(),frac),
                                  g.mul(inner.P_result(),g.sub(mpfr(1),g.div(mpfr(1),mpfr(frac))))+p_loss)
    if outer==None:
        return None
    return RSInnerOuter(inner.n,inner.k,inner.d,outer.n,outer.k,outer.d,cw.P_error(),p_loss)
    
def compare_rs_rate(L=[l for l in range(50,1000,20)],
                    pse=[1e-3], slope=[1e-3],
                    PL=lambda L, c, slope: g_pow(g.mul(mpfr(slope),mpfr(L)),mpfr(c))):

    XY = []
    XRE = []
    Label = []
    for s in slope:
        for p in pse:
            y = []
            re = []
            x = []
            Label.append( "RS p_break/nt={:1.1e} p_sym={:1.1e}".format(float(s),float(p)) )
            for l in L:
                rs = get_rs_at_res(1e9,l,mpfr("1e-14"),p,PL(l,1,s))
                if rs == None or rs.P_result() > mpfr("1e-14"):
                    continue
                x.append(l)
                y.append(rs.R_index(1e9))
                re.append(g.log10(rs.P_result()))
                
            XY.append( [x,y] )
            XRE.append( [x,re] )
    
    fig, ax = plt.subplots()
    plt.ylabel("Information Density")
    plt.xlabel("Strand Length (nt)")
    #ax.grid(True)
    for [x,y],label in zip(XY,Label):
        plt.plot(x,y,label=label)
    plt.legend()
    plt.show()
    dumpXY("compare_rs_rate",XY,Label)

    # fig, ax = plt.subplots()
    # plt.ylabel("Decoding Error Probability")
    # plt.xlabel("Strand Length (nt)")
    # #ax.grid(True)
    # for [x,y],label in zip(XRE,Label):
    #     plt.plot(x,y,label=label)
    # plt.legend()
    # plt.show()

    #dump(filename,X,"Length",Y,Label)

def compare_density_vs_loss(L=[l for l in range(50,1000,20)],
                    pse=[1e-3], slope=[1e-3],
                    PL=lambda L, c, slope: g_pow(g.mul(mpfr(slope),mpfr(L)),mpfr(c))):

    XY = []
    XRE = []
    Label = []
    for l in L:
        for p in pse:
            y = []
            re = []
            x = []
            Label.append( "RS L={} p_sym={:1.1e}".format(l,float(p)) )
            for s in slope:
                rs = get_rs_at_res(1e9,l,mpfr("1e-14"),p,PL(l,1,s))
                if rs == None or rs.P_result() > mpfr("1e-14"):
                    continue
                x.append(g.log10(s))
                y.append(rs.R_index(1e9))
                re.append(g.log10(rs.P_result()))
                
            XY.append( [x,y] )
            XRE.append( [x,re] )
    
    fig, ax = plt.subplots()
    plt.ylabel("Information Density")
    plt.xlabel("Log Probability of Strand Loss")
    #ax.grid(True)
    ax.set_xlim(-2.5,-8.5)
    for [x,y],label in zip(XY,Label):
        plt.plot(x,y,label=label)
    plt.legend()
    plt.show()
    dumpXY("compare_density_loss",XY,Label)

    # fig, ax = plt.subplots()
    # plt.ylabel("Decoder Error Probability")
    # plt.xlabel("Strand Length (nt)")
    # ax.set_xlim(-1,-8)
    # #ax.grid(True)
    # for [x,y],label in zip(XRE,Label):
    #     plt.plot(x,y,label=label)
    # plt.legend()
    # plt.show()
        

def rs_code_with_copy(n,d,copy,p_sym,p_erasure):
    loss = g_pow(p_erasure,copy)
    rs = RSCode(n,n-d+1,d,p_sym,loss)
    return rs
    
    
def compare_copy_to_code(p_err):
    label = []
    error = []
    rate = []
    rs = rs_code_with_copy(255,3,4,mpfr("0"),mpfr(p_err))
    
    label.append( rs.label() + " copies=4 " )
    error.append( g.log10(rs.P_result()) )
    rate.append( rs.R() / 4 )

    rs =get_outer_at_target(rs.P_result(),mpfr("0"),mpfr(p_err))
    label.append( rs.label() + " copies=1" )
    error.append( g.log10(rs.P_result()) )
    rate.append( rs.R() )

    for l,e,r in zip(label,error,rate):
        print ("{}: {} {}".format(l,e,r))
    
def make_plots():
    rs_sweep_length(n=255,
                k=[223],Length=[_ for _ in range(50,1025,25)],
                p_se=[mpfr(1e-2),mpfr(1e-3)],
                copies = [1,3,10],
                PL=lambda l,c: g_pow(g.mul(mpfr("0.0005"),mpfr(l)),mpfr(c)))

    compare_rs_rate(slope=[1e-3,5e-4,1e-4],L=[l for l in range(50,1020,20)])
    compare_density_vs_loss(slope=[ 10**(_/10.0) for _ in range(-80,-9, 10) ],L=[l for l in range(200,1100,200)])
    
#test:
#plot_binomials(256)    

if __name__ == "__main__":    
    make_plots()
