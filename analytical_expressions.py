import sympy as sy, numpy as np

# Symbols
C 				= sy.Symbol('C')
a 				= sy.Symbol('a')
alpha 			= sy.Symbol(r'\alpha')
beta 			= sy.Symbol(r'\beta')
omega			= sy.Symbol(r'\omega')
x1				= sy.Symbol('x_1')
y1				= sy.Symbol('y_1')
x2				= sy.Symbol('x_2')
y2				= sy.Symbol('y_2')
r12_abs			= sy.Symbol(r'|\mathbf{r}_1 - \mathbf{r}_2|')
r12_shorthand 	= sy.Symbol(r'r_{12}')
r1_squared 		= sy.Symbol(r'r_1^2')
r2_squared 		= sy.Symbol(r'r_2^2')
factor_of_half 	= sy.Symbol(r'\frac{1}{2}')

# Factors-to-be-differentiated list
diff_list = [x1,y1,x2,y2]

# Undifferentiated expressions
r12 			= sy.sqrt((x1-x2)**2 + (y1-y2)**2)
expression 		= C * sy.exp(-omega*alpha/2.*(x1**2 + y1**2 + x2**2 + y2**2) + a * r12 / (1 + beta * r12))
subexpression 	= -omega*alpha/2.*(x1**2 + y1**2 + x2**2 + y2**2) + a*r12/(1+beta*r12)
last_pt_expr	= a*r12/(1+beta*r12)
pt_2b_expr		= a*(x1-x2)/(r12*(1+beta*r12))
denom_r12beta	= 1/(r12*(1+beta*r12))

def full_derivative():
	# Differentiated expressions
	dpsi = [sy.diff(expression,x) for x in diff_list]
	ddpsi = [sy.diff(f,x) for f,x in zip(dpsi,diff_list)]
	# Expanding and adding expressions
	ddpsi_simply 	= sum(sy.expand(expr) for expr in ddpsi)
	ddpsi_simply	= ddpsi_simply / expression
	ddpsi_simply 	= sy.simplify(ddpsi_simply)
	# Substituting back into expression
	ddpsi_subs 		= ddpsi_simply.subs(r12, r12_shorthand)
	ddpsi_subs2		= ddpsi_subs.subs(x1**2 + y1**2, r1_squared)
	ddpsi_subs3		= ddpsi_subs2.subs(x2**2 + y2**2, r2_squared)
	return dpsi_subs3

# First derivate of 1/r_12
def derivate_r12(x):
	# Derivates one over r12
	dr12 			= sy.diff(1./r12,x)
	# Substituting for r12_shorthand
	dr12 			= dr12.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	return dr12

# First derivative of subexpression
def derivate_exponent_contents(x):
	ddx_exp			= sy.diff(subexpression,x)
	# Substituting back into expression
	ddx_exp			= ddx_exp.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	ddx_exp 		= ddx_exp.subs(x1**2+y1**2,r1_squared)
	ddx_exp 		= ddx_exp.subs(x2**2+y2**2,r1_squared)
	return ddx_exp

# Second derivative of subexpression
def double_derivate_exponent_contents(x):
	ddx_exp 		= sy.diff(subexpression,x)
	ddxddx_exp		= sy.diff(ddx_exp,x)
	# Substituting back into expression
	ddxddx_exp		= ddxddx_exp.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	ddxddx_exp 		= ddxddx_exp.subs(x1**2 + y1**2,r1_squared)
	ddxddx_exp 		= ddxddx_exp.subs(x2**2 + y2**2,r1_squared)
	return ddxddx_exp

def dx_psi(x):
	ddx_psi 		= sy.diff(expression,x)
	# ddx_psi 		= ddx_psi.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	return ddx_psi

def dxdx_psi(x):
	ddxddx_psi 		= sy.diff(expression,x)
	ddxddx_psi 		= sy.diff(ddxddx_psi,x)
	# ddxddx_psi 		= ddxddx_psi.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	return ddxddx_psi

def last_pt_derivate(x):
	dx_f 			= sy.diff(last_pt_expr, x)
	# dxdx_f 			= sy.diff(dx_f)
	dx_f 			= dx_f.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	return dx_f

def derivate_pt_2b(x):
	dx_pt_2b_expr 	= sy.diff(pt_2b_expr, x)
	dx_pt_2b_expr	= dx_pt_2b_expr.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	return dx_pt_2b_expr

def derivate_denom_r12beta(x):
	dx_return 		= sy.diff(denom_r12beta, x)
	dx_return 		= dx_return.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	return dx_return

def print_expressions(before,after,x):
	# Full expression:
	derivative 		= 'd^2'
	var 			= x
	print r'\frac{%s}{%s %s} \big[' % (derivative, derivative, var), sy.latex(before.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)), r'\big]'
	print r'=\\', sy.latex(after), r'\\'
	# print '\n\n\nRaw final: '
	# print after

var_to_derivate 	= x1

# final_expression 	= full_derivative()
# final_expression 	= derivate_r12(var_to_derivate)
# final_expression 	= derivate_exponent_contents(var_to_derivate)
# final_expression 	= double_derivate_exponent_contents(var_to_derivate)

final_expression 	= derivate_pt_2b(var_to_derivate)
print_expressions(pt_2b_expr, final_expression, var_to_derivate)

final_expression 	= derivate_denom_r12beta(var_to_derivate)
print_expressions(denom_r12beta, final_expression, var_to_derivate)

expr1 = r12*(1+beta*r12)
print_expressions(expr1, sy.diff(expr1,x1).subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand), x1)

expr2 = a*beta*(x1-x2)/(1+beta*r12)
print_expressions(expr2, sy.diff(expr2,x1).subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand), x1)

# final_expression	= dxdx_psi(var_to_derivate)
# final_expression	= final_expression.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
# final_expression	= final_expression.subs((x1 - x2)**2 + (y1 - y2)**2, r12_shorthand**2)
# final_expression	= final_expression.subs(x1**2 + y1**2,r1_squared)
# final_expression	= final_expression.subs(x2**2 + y2**2,r2_squared)
# final_expression	= final_expression/expression
# print_expressions(subexpression.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand), final_expression, var_to_derivate)

# final_expression 	= last_pt_derivate(var_to_derivate)
# print_expressions(last_pt_expr, final_expression,var_to_derivate)

# print_expressions(expression, final_expression, var_to_derivate)
