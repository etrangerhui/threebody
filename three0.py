import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve
import sys
ALWAYS=True
def getd(s, t, a):
    if ALWAYS:
        d = np.abs(s) * (t + a)
    else:
        d = np.abs(s) * t + a
    return d

def fgetline(fid):
    tline = fid.readline().strip()
    while tline == '' or tline == '\n' or tline == '\r':
        tline = fid.readline().strip()
    return tline

def readfile(filename):
    
    params =np.array( [], dtype=np.float64); 
    
    with open(filename, 'r') as f:
        try:
            tline = fgetline(f)
            parts = tline.split()
            nbody = int(parts[0])
            mcomp = int(parts[1])
            
            if not (2 <= nbody <= 30) or not (2 <= mcomp <= 3):
                raise ValueError("Invalid file format")
            
            
            a = np.zeros((nbody, mcomp), dtype=object)
            b = np.zeros((nbody, mcomp), dtype=object)
            c = np.zeros((nbody, mcomp), dtype=object)
            d = np.zeros((nbody, mcomp), dtype=object)
            
            tline = fgetline(f)
            parts = tline.split()
            a0 = np.array([float(parts[i]) for i in range(nbody)])
            m0 = np.ones(nbody)
            if len(parts) > nbody:
                m0 = np.array([float(parts[i]) for i in range(nbody, nbody * 2)])
            
            for i in range(nbody):
                for k in range(mcomp):
                    tline = fgetline(f)
                    parts = tline.split()
                    na = int(parts[0])
                    nd = int(parts[1])
                    a[i][k] = np.zeros(na, dtype=int)
                    b[i][k] = np.zeros(na, dtype=np.float64)
                    c[i][k] = np.zeros(nd, dtype=np.float64)
                    d[i][k] = np.zeros(nd, dtype=int)
                    
                    if na > 0:
                        tline = fgetline(f)
                        parts = tline.split()
                        a[i][k] = np.array([int(float(parts[m])) for m in range(na)])
                        if len(parts) > na:
                            for m in range(na, len(parts)):
                                    b[i][k][m-na] = np.float64(parts[m])
                    if nd > 0:
                        tline = fgetline(f)
                        parts = tline.split()
                        d[i][k] = np.array([int(float(parts[m])) for m in range(nd)])
                        if len(parts) > nd:
                            for m in range(nd, len(parts)):
                                    c[i][k][m-nd] = np.float64(parts[m])
                    params=np.concatenate((params,b[i][k][a[i][k]>0],c[i][k][d[i][k]>0]))
        
        except EOFError:
            raise IOError("Failed to read the file")
        except Exception as e:
            raise IOError(f"Error reading the file: {e}")    
    return a,b,c, d, a0, m0, params

def P(t, a, b,c,d, n, m, a0):
    p = 0.0    
    if len(a[n][m])> 0:
        p += np.sum(b[n][m] * np.sin(getd(a[n][m], t, a0[n]))) 
    if len(d[n][m]) > 0:
        p += np.sum(c[n][m] * np.cos(getd(d[n][m], t, a0[n])))    
    return p
def Q(t, a, b,c,d, n, m, a0):
    p = 0.0  
    if len(a[n][m]) > 0:
        p += np.sum(a[n][m]*b[n][m] * np.cos(getd(a[n][m], t, a0[n]))) 
    if len(d[n][m]) > 0:
        p -= np.sum(d[n][m]*c[n][m] * np.sin(getd(d[n][m], t, a0[n])))    
    return p
def R(t, a, b,c,d, n, m, a0):
    p = 0.0  
    if len(a[n][m]) > 0:
        p -= np.sum(a[n][m]*a[n][m]*b[n][m] * np.sin(getd(a[n][m], t, a0[n]))) 
    if len(d[n][m]) > 0:
        p -= np.sum(d[n][m]*d[n][m]*c[n][m] * np.cos(getd(d[n][m], t, a0[n])))    
    return p

def myFun(t, a,b,c, d, a0, n, m, m0,v,j):
    f = 0.0
    p, q = len(a), len(a[0])
    for nn in range(p):
        if n != nn:
            r2 = 0.0
            x = np.zeros(q)
            for i in range(q):
                u = (m + i) % q
                x[i] = P(t,a,b,c,d, nn, u, a0) - P(t, a,b,c, d, n, u, a0)
                r2 += x[i]**2
            
            r3 = r2 * np.sqrt(r2)
            if v<0 and j<0:
                f += m0[nn] * x[0] / r3
            else:
                r5 = r2 * r3
                j_mod = (j - m + q) % q
                
                if j_mod == 0:
                    if nn == v:
                        f += m0[nn] * (r2 - 3 * x[0] * x[0]) / r5
                    elif n == v:
                        f -= m0[nn] * (r2 - 3 * x[0] * x[0]) / r5
                elif j_mod == 1:
                    if nn == v:
                        f -= m0[nn] * 3 * x[0] * x[1] / r5
                    elif n == v:
                        f += m0[nn] * 3 * x[1] * x[0] / r5
                elif j_mod == 2:
                    if nn == v:
                        f -= m0[nn] * 3 * x[0] * x[2] / r5
                    elif n == v:
                        f += m0[nn] * 3 * x[2] * x[0] / r5                
    return f


def update(params,a, b, c,d):
    m, n = a.shape[0], a.shape[1]  
    u = 0  
    for i in range(m):
        for j in range(n):
            for k in range(len(a[i][j])):
                if(a[i][j][k]>0):
                   b[i][j][k] = params[u]
                   u += 1  
            for k in range(len(d[i][j])):
                if(d[i][j][k]>0):
                   c[i][j][k] = params[u]
                   u += 1 
    return b, c

def func(params, a,b, c, d, a0, m0, PI):
    m, n = a.shape[0], a.shape[1]  
    b, c = update(params, a,b, c, d,)  
    f = []  
    u, v = -1, -1
    for i in range(m):
        for j in range(n):
            k = len(a[i][j])  
            for p in range(k):
                s = a[i][j][p];
                if s<=0:continue;
                w, _ = quad(lambda t: myFun(t, a,b,c, d, a0, i, j, m0,u,v) * np.sin(getd(s, t, a0[i])), 0, 2 * PI)
                f.append(-s**2 * b[i][j][p] - w / PI)
            k = len(d[i][j])  
            for p in range(k):
                s = d[i][j][p]
                if s<=0:continue;
                w, _ = quad(lambda t: myFun(t, a,b,c, d, a0, i, j, m0,u,v) * np.cos(getd(s, t, a0[i])), 0, 2 * PI)
                f.append(-s**2 * c[i][j][p]- w / PI)
    s=0;
    k=int(2*np.pi/0.1)
    for p in range(k):
        t=0.1+p*0.1
        for i in range(m):
            for j in range(n):
                ds=myFun(t, a,b,c, d, a0, i, j, m0,u,v)-R(t, a, b,c,d, i, j, a0)
                s+=ds*ds       
    return f,s**0.5

def Dfunc(params, a,b, c, d, a0, m0, PI):
    m, n = a.shape[0], a.shape[1]  
    b, c = update(params, a,b,c, d)     
    f = []  

    for i in range(m):
        for j in range(n):            
            for k in range(len(a[i][j])):
                f0 = []
                s = a[i][j][k] 
                if s<=0:continue;
                for u in range(m):
                    for v in range(n):                        
                        for p in range(len(a[u][v])):                            
                            q = a[u][v][p]  
                            if q<=0:continue;
                            w, _ = quad(lambda t: myFun(t, a,b,c, d, a0, i, j, m0, u, v) * np.sin(getd(s, t, a0[i])) * np.sin(getd(q, t, a0[u])), 0, 2 * PI)
                            if i == u and j == v and p == k:
                                f0.append(-s**2 - w/PI)
                            else:
                                f0.append(-w/PI)                        
                        for p in range(len(d[u][v])):
                            q = d[u][v][p]
                            if q<=0:continue;
                            w, _ = quad(lambda t: myFun(t, a,b,c, d, a0, i, j, m0, u, v) * np.sin(getd(s, t, a0[i])) * np.cos(getd(q, t, a0[u])), 0, 2 * PI)
                            f0.append(-w/PI)
                f.append(f0)
            
           
            for k in range(len(d[i][j])):
                f0 = []
                s = d[i][j][k] 
                if s<=0:continue;
                for u in range(m):
                    for v in range(n):
                        for p in range(len(a[u][v])):
                            q = a[u][v][p] 
                            if q<=0:continue;
                            w, _ = quad(lambda t: myFun(t, a,b,c, d, a0, i, j, m0, u, v) * np.cos(getd(s, t, a0[i])) * np.sin(getd(q, t, a0[u])), 0, 2 * PI)
                            f0.append(-w/PI)
                        for p in range(len(d[u][v])):
                            q = d[u][v][p]
                            if q<=0:continue;
                            w, _ = quad(lambda t: myFun(t, a,b,c, d, a0, i, j, m0, u, v) * np.cos(getd(s, t, a0[i])) * np.cos(getd(q, t, a0[u])), 0, 2 * PI)
                            if i == u and j == v and p == k:
                                f0.append(-s**2 - w/PI)
                            else:
                                f0.append(-w/PI)
                f.append(f0)
    return f



def solver(x0, a,b, c, d, a0, m0, PI, tol, max_iter):
    x = x0
    g,diff = func(x, a,b, c, d, a0, m0, PI)
    H = Dfunc(x, a,b, c, d, a0, m0, PI)
    iter = 0
    while iter < max_iter:
        if np.linalg.norm(g) < tol:
            root = x
            print(f'Converged to root after {iter} iterations.')
            return root
        p = solve(H, g)
        x_new = x - p
        g_new,diff = func(x_new, a,b, c, d, a0, m0, PI)
        H = Dfunc(x_new, a,b, c, d, a0, m0, PI)
        x = x_new
        g = g_new
        print(np.linalg.norm(p), np.linalg.norm(g),diff)
        iter += 1
    root = x
    print('Maximum iterations reached. Solution may not be accurate.')
    return root

if __name__=="__main__":
    G = 6.67430e-11
    if len(sys.argv)>1:
        a, b,c,d, a0, m0, params=readfile(sys.argv[1]);
    else:
        a, b,c,d, a0, m0, params=readfile('cross.txt');
    if len(sys.argv)>2 and sys.argv[2].lower() in ('yes', 'true', 't', '1'):ALWAYS = True
    if len(sys.argv)>2 and sys.argv[2].lower() in ('no', 'false', 'f', '0'):ALWAYS = False
    params = solver(params, a, b,c,d, a0, m0, np.pi, 1.e-13, 100)

    with open('test.unv', 'w') as fp:
        fp.write(str(len(m0)) +' 1000 1.e-12 \n')
        t = 0
        b, c = update(params, a,b,c,d)
        for i in range(len(m0)):
            p1 = P(t, a, b,c,d, i, 0, a0)
            p2 = Q(t, a, b,c,d, i, 0, a0)
            p3 = P(t, a, b,c,d, i, 1, a0)
            p4 = Q(t, a, b,c,d, i, 1, a0)
        
 
            fp.write(f'{str(p1)} {str(p3)} 0.0 {str(m0[i]/G)} 0.01 {str(p2)} {str(p4)} 0.0 {i+1}\n')
        nbody=len(a);mcomp=len(a[0]);
        for i in range(nbody):
            for j in range(mcomp):
                for k in range(len(a[i][j])):
                    if ALWAYS:
                        fp.write(f"{str(b[i][j][k])}*sin({str(abs(a[i][j][k]))}*(t+{str(a0[i])}))+")
                    else:
                        fp.write(f"{str(b[i][j][k])}*sin({str(abs(a[i][j][k]))}*t+{str(a0[i])})+")
                for k in range(len(d[i][j])):
                    if ALWAYS:
                        fp.write(f"{str(c[i][j][k])}*cos({str(abs(d[i][j][k]))}*(t+{str(a0[i])}))+")
                    else:
                        fp.write(f"{str(c[i][j][k])}*cos({str(abs(d[i][j][k]))}*t+{str(a0[i])})+")
                fp.write(f"[{i+1}][{j+1}]\n")
        fp.write(f"{nbody} {mcomp}\n")
        for i in range(nbody):
            fp.write(f"{str(a0[i])} ")
        for i in range(nbody):
            fp.write(f"{str(m0[i])} ")
        fp.write("\n")
        for i in range(nbody):
            for j in range(mcomp):
                fp.write(f"{len(a[i][j])} {len(d[i][j])}   #[{i+1},{j+1}]\n")
                for k in range(len(a[i][j])):
                    fp.write(f"{str(a[i][j][k])} ")
                for k in range(len(b[i][j])):
                    fp.write(f"{str(b[i][j][k])} ")
                fp.write("\n")
                for k in range(len(d[i][j])):
                    fp.write(f"{str(d[i][j][k])} ")
                for k in range(len(c[i][j])):
                    fp.write(f"{str(c[i][j][k])} ")
                fp.write("\n")
