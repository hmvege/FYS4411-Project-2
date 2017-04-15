import sympy as sy, numpy as np

def return_sum(expr_list):
	sym_obj = 0
	for sym in expr_list:
		sym_obj += sym
	return sym_obj

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

# Factors-to-be-differentiated list
diff_list = [x1,y1,x2,y2]

# Undifferentiated expressions
r12 			= sy.sqrt((x1-x2)**2 + (y1-y2)**2)
expression 		= C * sy.exp(-omega/2.*(x1**2 + y1**2 + x2**2 + y2**2) + a * r12 / (1 + beta * r12))
subexpression 	= -omega/2.*(x1**2 + y1**2 + x2**2 + y2**2) + a*r12/(1+beta*r12)

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
def derivate_r12():
	# Derivates one over r12
	dr12 = sy.diff(1./r12,x1)
	# Substituting for r12_shorthand
	dr12 = dr12.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	return dr12

# First derivative of subexpression
def derivate_exponent_contents():
	ddx_psi			= sy.diff(subexpression,x1)
	# Substituting back into expression
	ddx_psi			= ddx_psi.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	return ddx_psi

# Second derivative of subexpression
def double_derivate_exponent_contents():
	ddx_psi 			= derivate_exponent_contents()
	ddxddx_psi			= sy.diff(ddx_psi,x1)
	# Substituting back into expression
	ddxddx_psi			= ddxddx_psi.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)
	return ddx_psi



final_expression = derivate_exponent_contents()
final_expression = double_derivate_exponent_contents()

def dx_exponent():
	ddx_psi = sy.diff(expression,)

def ddx_exponent():


# final_expression = full_derivative()
# final_expression = derivate_r12()

# ddx_subs 		= ddx_simpl.subs(r12,r12_shorthand)
# ddx_subs2		= ddx_subs.subs(x1**2 + y1**2, r1_squared)
# ddx_subs3 	= ddx_subs2.subs(x2**2 + y2**2, r2_squared)
# final_expression = ddx_subs3

# final_expression = sy.simplify(ddpsi_subs3)

# Full expression:
print 'Raw final: '
print final_expression
print '\n\n\nLatex final: '
print sy.latex(final_expression)
