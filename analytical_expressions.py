import sympy as sy, numpy as np

class DerFinder:
	def __init__(self, expression):
		self.expression = expression
		self.substitutions = []

	def add_subs(self, subs_expr):
		"""
		Arguments
			subs_expr: 	Takes a list where each cell is a list of length 2
		"""
		self.substitutions = subs_expr

	def _derivate(self, vars_to_derivate):
		"""
		Internal function used for derivating. Takes list of derivatives and adds them together.
		"""
		self.expression_differentiated = sum(sy.diff(self.expression, x) for x in vars_to_derivate)

	def _substitute(self, expr):
		"""
		Internal function used when substituting expression
		"""
		if len(self.substitutions) == 0: return expr
		for expr_element in self.substitutions:
			expr = expr.subs(expr_element[0], expr_element[1])
		return expr

	def __call__(self, vars_to_derivate, print_expr=True, simplify=False):
		"""
		When called, find the derivative, substitutes and prints the expression on latex form
		"""
		assert type(vars_to_derivate) is list, "Derivations must be contained in a list: %s is not valid" % vars_to_derivate
		self._derivate(vars_to_derivate)
		self.expression_differentiated = self._substitute(self.expression_differentiated)
		self.expression = self._substitute(self.expression)
		if simplify: self.expression_differentiated = sy.simplify(self.expression_differentiated)
		if print_expr:
			if len(vars_to_derivate) == 1:
				derivative = 'd'
			else:
				derivative = 'd^%d' % (len(vars_to_derivate))
			print r'\frac{%s}{%s %s} \left(' % (derivative, derivative, sum(i for i in vars_to_derivate)), sy.latex(self.expression), r'\right)'
			print r'=\\', sy.latex(self.expression_differentiated), r'\\'
		return self.expression_differentiated

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
r12r12		 	= sy.Symbol(r'r_{12}^2')
r1_squared 		= sy.Symbol(r'r_1^2')
r2_squared 		= sy.Symbol(r'r_2^2')
factor_of_half 	= sy.Symbol(r'\frac{1}{2}')

# Substitution lists
substitutions_list = [	[sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand],
						[x1**2 + y1**2, r1_squared],
						[x2**2 + y2**2, r2_squared],
						[(x1-x2)**2 + (y1-y2)**2, r12r12]]

# Factors-to-be-differentiated list
diff_list = [x1,y1,x2,y2]

# Undifferentiated expressions
r12 			= sy.sqrt((x1-x2)**2 + (y1-y2)**2)
expression 		= C * sy.exp(-omega*alpha/2.*(x1**2 + y1**2 + x2**2 + y2**2) + a * r12 / (1 + beta * r12))
exponent 		= -omega*alpha/2.*(x1**2 + y1**2 + x2**2 + y2**2) + a*r12/(1+beta*r12)
last_pt_expr	= a*r12/(1+beta*r12)
pt_2b_expr		= a*(x1-x2)/(r12*(1+beta*r12))
denom_r12beta	= 1/(r12*(1+beta*r12))

def derivate_1_over_r12(der):
	# First derivate of 1/r_12
	dr12 = DerFinder(1./r12)
	dr12.add_subs([[sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand]])
	dr12(der)

def derivate_exponent(der):
	d_exp = DerFinder(exponent)
	d_exp.add_subs(substitutions_list)
	expr = d_exp(der,simplify=True)
	# print sy.latex(sy.simplify(expr*(-omega*alpha*x1) - omega*alpha))

def derivate_beta_frac(der):
	d_beta = DerFinder(1/(1+beta*r12))
	d_beta.add_subs(substitutions_list)
	d_beta(der,simplify=True)

def derivate_frac(der):
	d_frac = DerFinder(a*(x1-x2)/(r12*(1+beta*r12)**2))
	d_frac.add_subs(substitutions_list)
	d_frac(der,simplify=False)

def double_derivate_exponent(der):
	d_exp = DerFinder(exponent)
	# d_exp.add_subs(substitutions_list)
	d_expr = d_exp(der,simplify=False,print_expr=True)
	dd_expr = DerFinder(d_expr)
	dd_expr.add_subs(substitutions_list)
	dd_expr(der)

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

def dpsi_beta(der):
	expr = DerFinder(expression)
	expr.add_subs(substitutions_list)
	dPsiBeta = expr([der],simplify=True,print_expr=True)


def print_expressions(before,after,x):
	# Full expression:
	derivative 		= 'd^2'
	var 			= x
	print r'\frac{%s}{%s %s} \left[' % (derivative, derivative, var), sy.latex(before.subs(sy.sqrt((x1-x2)**2 + (y1-y2)**2), r12_shorthand)), r'\big]'
	print r'=\\', sy.latex(after), r'\\'
	# print '\n\n\nRaw final: '
	# print after

vars_to_derivate = [x2]

def main():
	# derivate_1_over_r12()
	# derivate_exponent(vars_to_derivate)
	# derivate_beta_frac(vars_to_derivate)
	# derivate_frac(vars_to_derivate)
	# double_derivate_exponent(vars_to_derivate)
	dpsi_beta(alpha)

if __name__ == '__main__':
	main()