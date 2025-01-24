from mpmath import mp
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import deque
import copy
import re

mp.dps=50
#    a2 =mp.mpf( 0.25);
#    a3 =mp.mpf( 0.375);
#    a4 = mp.mpf(12)/mp.mpf(13);
#    a6 =mp.mpf( 0.5);
half=mp.mpf( 0.5      )     
G = mp.mpf(6.67430e-11)
b21 =mp.mpf( 0.25      )            
b31 =mp.mpf( 3)/mp.mpf( 32      )   
b32 =mp.mpf( 9)/mp.mpf( 32      )   
b41 =mp.mpf( 1932)/mp.mpf( 2197 )   
b42 =mp.mpf( -7200)/mp.mpf( 2197)   
b43 =mp.mpf( 7296)/mp.mpf( 2197 )   
b51 =mp.mpf( 439)/mp.mpf( 216   )   
b52 =mp.mpf( -8        )   
b53 =mp.mpf( 3680)/mp.mpf( 513  )   
b54 =mp.mpf( -845)/mp.mpf( 4104 )   
b61 =mp.mpf( -8)/mp.mpf( 27     )   
b62 =mp.mpf( 2         )  ;
b63 =mp.mpf( -3544)/mp.mpf( 2565)   
b64 =mp.mpf( 1859)/mp.mpf( 4104 )   
b65 =mp.mpf( -11)/mp.mpf( 40    )   
c1 = mp.mpf(25)/mp.mpf( 216     )   
c3 = mp.mpf(1408)/mp.mpf( 2565  )   
c4 = mp.mpf(2197)/mp.mpf( 4104  )   
c5 = mp.mpf(-0.20      )     
d1 = mp.mpf(1)/mp.mpf( 360      )   
d3 = mp.mpf(-128)/mp.mpf( 4275  )   
d4 = mp.mpf(-2197)/mp.mpf( 75240)   
d5 = mp.mpf(0.02       )       
d6 = mp.mpf(2)/mp.mpf( 55       )  
def apply_op(a, b, op) :
    if op == '+':
        return a + b
    elif op == '-':
        return a - b
    elif op == '*':
        return a * b
    elif op == '/':
        if b == mp.mpf(0):
            raise ValueError("Division by zero")
        return a / b
    elif op == '^':
        return a ** b
    else:
        raise ValueError(f"Unknown operator: {op}")

def evaluate(expression) :
    values = []  # 用于存储数值的栈
    ops = []     # 用于存储操作符的栈

    # 处理表达式字符串，替换特殊字符并添加空格以便于分词
    expression = re.sub(r'^-', '0-', expression)
    expression = re.sub(r'\+', ' + ', expression)
    expression = re.sub(r'e \+ ', 'e+', expression)
    expression = re.sub(r'E \+ ', 'e+', expression)
    expression = re.sub(r'-', ' - ', expression)
    expression = re.sub(r'e - ', 'e-', expression)
    expression = re.sub(r'E - ', 'e-', expression)
    expression = re.sub(r'\*', ' * ', expression)
    expression = re.sub(r'/', ' / ', expression)
    expression = re.sub(r'\(', ' ( ', expression)
    expression = re.sub(r'\)', ' ) ', expression)
    expression = re.sub(r'\^', ' ^ ', expression)

    # 分词
    tokens = re.findall(r'\S+', expression)

    for token in tokens:
        if token[0].isdigit() and len(token)>=1:
            # 如果是数字（包括负数），则转换为mp.mpf并压入栈中
            values.append(mp.mpf(token))
        elif token == '(':
            # 如果是左括号，则压入栈中
            ops.append(token)
        elif token == ')':
            # 如果是右括号，则计算直到遇到左括号
            while ops and ops[-1] != '(':
                b = values.pop()
                a = values.pop()
                op = ops.pop()
                values.append(apply_op(a, b, op))
            # 弹出左括号
            if ops:
                ops.pop()
        else:
            # 如果是操作符，则根据优先级处理
            while ops and precedence(ops[-1]) >= precedence(token):
                b = values.pop()
                a = values.pop()
                op = ops.pop()
                values.append(apply_op(a, b, op))
            ops.append(token)

    # 处理栈中剩余的操作符
    while ops:
        b = values.pop()
        a = values.pop()
        op = ops.pop()
        values.append(apply_op(a, b, op))

    # 表达式应该只有一个结果
    if len(values) != 1:
        raise ValueError("Invalid expression")

    return values[0]

# 定义操作符的优先级
def precedence(op) :
    if op in ('+', '-'):
        return 1
    if op in ('*', '/'):
        return 2
    if op == '^':
        return 3
    return 0





# 物体类
class Object:
    def __init__(self, x, y,vx,vy, mass, radius):
        self.position = [x, y]
        self.velocity = [vx, vy]
        self.mass = mass
        self.radius = radius

    def __mul__(self, other):
        return Object(self.position[0] * other,self.position[1] * other,self.velocity[0] * other,self.velocity[1] * other,self.mass,self.radius)
        
    def __add__(self, other):   
        return Object(self.position[0] + other.position[0],self.position[1]+ other.position[1],self.velocity[0] + other.velocity[0],self.velocity[1]+ other.velocity[1],self.mass,self.radius)
          

# 计算两个物体之间的引力
def calculate_gravity(obj1, obj2):
    dx=[obj1.position[0] - obj2.position[0],obj1.position[1] - obj2.position[1]];
    distance = dx[0]*dx[0]+dx[1]*dx[1]
    if distance == 0:
        return [0, 0] # 避免除以零
    force_magnitude = G * obj1.mass * obj2.mass / distance    
    return [dx[0]*force_magnitude  / distance**half,dx[1]*force_magnitude  / distance**half]
def calculate_acc(obj2, objs):
    a=[0,0]
    for i in range(len(objs)):
        obj1 = objs[i]
        dx=[obj1.position[0] - obj2.position[0],obj1.position[1] - obj2.position[1]];
        distance = dx[0]*dx[0]+dx[1]*dx[1]
        if distance == 0:continue
        force_magnitude = G * obj1.mass  / distance

        a= [a[0]+dx[0]*force_magnitude  / distance**half,a[1]+dx[1]*force_magnitude  / distance**half]
    
    return a
def fgetline(fid):
    tline = fid.readline().strip()
    while tline == '' or tline == '\n' or tline == '\r':
        tline = fid.readline().strip()
    return tline

def update_objects(objects, dt):
    for i in range(len(objects)):
        for j in range(len(objects)):
            if i==j:continue
            obj1 = objects[i]
            obj2 = objects[j]
            
            # 计算引力
            force = calculate_gravity(obj1, obj2)
            acceleration = force / obj1.mass
            obj1.velocity += acceleration * dt
            
    for i in range(len(objects)):
        obj1 = objects[i]
        obj1.position += obj1.velocity * dt
        tracex[i].append(obj1.position[0])  
        tracey[i].append(obj1.position[1])

def getdy(obj):
    objects=[]
    for i in range(len(obj)):
       obj1 = obj[i]
       objects+=[copy.deepcopy(obj1)]
       a=calculate_acc(obj1, obj)   
       objects[i].position=obj1.velocity
       objects[i].velocity=a
    return objects
def norm(obj):
    objects=mp.mpf(0)
    for i in range(len(obj)):
       obj1 = obj[i]
       for j in range(len(obj1.position)):
           objects+= obj1.position[j]*obj1.position[j]
       for j in range(len(obj1.velocity)):
           objects+= obj1.velocity[j]*obj1.velocity[j]
    return objects**half
def getscale(objs,ratio):
    obj = copy.deepcopy(objs[0])
    if ratio[0]!=1:
        for j in range(len(objs[0])): 
            obj[j]*=ratio[0]
    for i in range(1,len(objs)):
        for j in range(len(objs[i])):
            obj[j]+=objs[i][j]*ratio[i]
    return obj
def  Runge_Kutta(y,h):
 
    
  #  h2 = a2 * h; h3 = a3 * h; h4 = a4 * h; h6 = a6 * h;
    
    k1 = getdy(y);                                                               #x0
    k2 = getdy(getscale([y,k1],[1, h * b21 ]));                                  #x0+h2
    k3 = getdy(getscale([y,k1,k2] ,[1,h *b31,h*b32]) );                            #x0+h3
    k4 = getdy(getscale([y,k1,k2,k3] ,[1, h *b41,h*b42, h*b43]));#                   #x0+h4 
    k5 = getdy(getscale([y,k1,k2,k3,k4] , [1,h *b51, h*b52, h*b53,h*b54]));             #x0+h 
    k6 = getdy(getscale([y,k1,k2,k3,k4,k5] , [1,h*b61,h*b62, h*b63,h*b64,h*b65]));      #x0+h6
    y = getscale([y,k1,k3,k4,k5],  [1,h * c1,h * c3,h * c4,h * c5]);
    out = getscale([k1,k3,k4,k5,k6],[d1,d3,d4,d5,d6]);
    return y,norm(out)

def ODE_RK4_hyh2(y0,h):
    y1,err=Runge_Kutta(y0,h);
    if err>mp.mpf(1.e-12):
        y1=ODE_RK4_hyh2(y0,h*half);
        y1=ODE_RK4_hyh2(y1,h*half);   
    return y1


objects=[]
with open('test1.unv', 'r') as f:
    try:
        tline = fgetline(f)
        parts = tline.split()
        nbody = int(parts[0])
        for i in range(nbody):
            tline = fgetline(f)
            parts = tline.split()
            objects+=[Object(evaluate(parts[0]), evaluate(parts[1]),evaluate(parts[5]), evaluate(parts[6]), evaluate(parts[3]), float(parts[4])) ]
    except EOFError:
        raise IOError("Failed to read the file")
    except Exception as e:
        raise IOError(f"Error reading the file: {e}")
dt = mp.mpf(2*3.1415926535897/100)  
num_frames = 1
tracex=[];tracey=[];
for i in range(nbody):
    tracex+= [deque(maxlen=500)]
    tracey+= [deque(maxlen=500)]

lines = []
fig, ax = plt.subplots()

ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')
scatter = ax.scatter(range(nbody), range(nbody),c=range(nbody))  # 初始化为空散点图
sizes = [o.radius*10000 for o in objects]
scatter.set_sizes(sizes)
for _ in range(nbody):
    line, = ax.plot([], [])  # 创建线条对象，但不立即绘制
    lines.append(line)  # 将线条添加到列表中




def animate(frame):
    global objects
    for _ in range(num_frames): 
        objects=ODE_RK4_hyh2(objects, dt)
        for i in range(len(objects)):
            tracex[i].append(objects[i].position[0])  
            tracey[i].append(objects[i].position[1])
    data=[];        
    for i in range(len(objects)):
        data.append(objects[i].position) 
    scatter.set_offsets(data)
    
    for  i,line in enumerate(lines):
        line.set_data(tracex[i], tracey[i]) 
    
    return [scatter]+lines

# 创建动画
ani = FuncAnimation(fig, animate, frames=num_frames, interval=100, blit=True)


plt.show()