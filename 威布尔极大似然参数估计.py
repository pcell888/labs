from sympy import *
import mpmath
mpmath.mp.dps = 15

x = symbols("x")

data = [153514.6227,187844.2601,247900.6145,275520.2084,282634.6126,316168.3311,316215.2135,317094.2585,366647.8734,414861.6434,414861.6434,414861.6434,414861.6434]

n, r, T = 13, 9, data[-1]

expr1 = 0
for t in data[:r]:
    expr1 = expr1 + t**x*log(t)

expr2 = (n-r)*T**x * log(T)

expr3 = 0
for t in data[:r]:
    expr3 = expr3 + t**x

expr3 = expr3 + (n-r)*T**x

expr4 = 1/x

expr5 = 0
for t in data[:r]:
    expr5 = expr5 + log(t)

expr5 = (1/r)*expr5

expr = simplify( (((expr1+expr2)/expr3) - expr4) - expr5)

# 计算β结果
beta = nsolve(Eq(expr, 0), x, 5, verify=False)
print('β：', beta) # 3.12245362959167

# 计算η结果
yita = ( (1/r)* (sum(data[:r]) + (n-r)*T**beta) ) **(1/beta)
print('η：', yita) # 319972.733682525