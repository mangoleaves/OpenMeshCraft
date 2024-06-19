#pragma once

#include <cfloat>
#include <cmath>
#include <cstdint>

#include <algorithm>
#include <deque>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

typedef double   fpnumber;
constexpr double FPN_EPSILON = DBL_EPSILON;

// In bytes, includes all: input parameters, doubles, ints, ...
#define MAX_STACK_SIZE 16384
// This is the maximum in any case, but might be reduced to avoid stack
// overflows
#define MAX_STATIC_SIZE 128
// If the polynomial degree is higher than this, avoid FP arithmetic at all
#define MAX_WORTH_DEGREE 20

// If defined, exact functions check for underflows and try to fix
#define UNDERFLOW_GUARDING

// floating-point number type in generated code
extern std::string FT;
// interval number type (IT) in generated code
extern std::string IT;
// exact number type (ET) in generated code
extern std::string ET;
// boolean type in generated code
extern std::string Boolean;
// sign type in generated code
extern std::string Sign;
// API symbol in generated code
extern std::string API;
// point arrangement used by filter in generated code
extern std::string PntArr;
// Exit command in generated code
extern std::string Exit;

/// @brief print help information and exit.
void help_info();

/// @brief print error and exit.
void error(const char *s, int line);

/// @brief separate string line to some tokens with specified separator.
void tokenize(const std::string &line, std::vector<std::string> &tokens,
              char separator);

/// unit in the last place.
/// ulp(x) is the gap between the two finite floating-point numbers nearest x,
/// even if x is one of them (But ulp(NaN) is NaN).
fpnumber ulp(fpnumber d);

class Variable
{
public:
	// Properties ###############################################################

	std::string name; // Variable name

	// TRUE if this var is the part of an explicit point
	// example: this var is part of explicitPoint(p:px,py,pz)
	bool is_part_of_explicit;

	// TRUE if this var is the part of an implicit point
	// example: this var is part of implicitPoint(p:bx,by,bz,lx,ly,lz,ld)
	bool is_part_of_implicit;

	bool is_lambda;   // TRUE if name starts with "lambda"
	bool is_lambda_d; // TRUE if name starts with "lambda_d"
	bool is_beta;     // TRUE if name starts with "beta"

	bool is_used; // TRUE if this var is the operand of another var

	// Operands and operator ####################################################

	Variable *op1, *op2; // Operands (pointers to all_vars. nullptr if input val)
	char      op;        // Operation (=, +, -, *)

	// Error ####################################################################

	bool error_evaluated; // TRUE if error is evaluated in error propagation
	bool is_a_max;        // TRUE if magnitude is relevant for error bound

	int      size;         // Variable size (expansion max length)
	int      error_degree; // Forward error analysis: degree
	// fp: floating-point double
	fpnumber fp_error_bound; // Forward error analysis: FP bound on error
	fpnumber fp_value_bound; // Forward error analysis: FP bound on magnitude

	std::string actual_length; // Variable length (expansion actual length)

public:
	// Auxiliary type distinguisher
	struct ExplicitT
	{
	};
	struct ImplicitT
	{
	};

public: /* Constructors ******************************************************/
	// Plain declaration
	Variable(const std::string &s);

	// Plain explicit declaration
	Variable(std::string &s, ExplicitT);

	// Plain implicit declaration
	Variable(std::string &s, ImplicitT);

	// Plain assignment
	Variable(std::string &s, Variable *o1);

	// Binary operation
	Variable(std::string &s, Variable *o1, Variable *o2, char o);

public: /* Error definitions *************************************************/
	void clearError();

	void setError(int _error_degree, int _size, fpnumber _fp_error_bound,
	              fpnumber _fp_value_bound);

	void propagateError();

public: /* Checks ************************************************************/
	bool isInput() const { return op1 == nullptr && op2 == nullptr; }

	bool isParameter() const
	{
		return isInput() && name != "1" && name != "2" && !is_part_of_implicit &&
		       !is_part_of_explicit;
	}
};

/**
 * @brief LambdaVariable comes from declaration of implicit point.
 */
class LambdaVariable
{
public:
	// Type of the implicit point (e.g. implicitPoint2D/3D)
	std::string             point_type;
	// Name of the implicit point (e.g. p1)
	std::string             point_name;
	// Indices in all_vars of the results of getLambda() (e.g. l1x, l1y, d)
	std::vector<Variable *> output_pars;
	// dimension of the point
	int                     dim;

public:
	LambdaVariable(std::string &n);

	std::string get_type_string() const;

	std::string print_filtered() const;
	std::string print_interval() const;
	std::string print_expansion() const;
	std::string print_exact() const;
};

/**
 * @brief ExplicitVariable comes from declaration of explicit point.
 */
class ExplicitVariable
{
public:
	// Type of the explicit point (e.g., explicitPoint2D/3D)
	std::string             point_type;
	// Name of the explicit point (e.g. p)
	std::string             point_name;
	// coordinates
	std::vector<Variable *> coords;
	// dimension of the point
	int                     dim;

public:
	ExplicitVariable(std::string &n);

	std::string get_type_string() const;

	std::string get_vars() const;
	std::string get_coords() const;
};

class ErrorDefinition
{
public:
	// Name of the implicit point type (e.g., SSI, LPI, TPI)
	std::string name;
	// Short name of the implicit point type (e.g., S, L, T)
	std::string short_name;
	// dimension
	size_t      dim;
	// error definitions
	struct ErrorDef
	{
		int      size;         // Variable size (expansion max length)
		int      error_degree; // Forward error analysis: degree
		// fp: floating-point double
		fpnumber fp_value_bound; // Forward error analysis: FP bound on magnitude
		fpnumber fp_error_bound; // Forward error analysis: FP bound on error
	};
	std::vector<ErrorDef> pars; // [x,y,z,d] in 3D, or [x,y,d] in 2D

public: /* Constructors ******************************************************/
	void parseOneErrorDef(ErrorDef &e, std::string &tokens);
};

class Predicate
{
public:
	bool append;           // TRUE to append code to an existing file
	bool output_filtered;  // TRUE to output filtered function
	bool output_interval;  // TRUE to output interval function
	bool output_exact;     // TRUE to output exact function
	bool output_expansion; // TRUE to output expansion function

	std::string support_operators = "+-*";

	// List of all the variables
	std::deque<Variable>   all_vars;
	// List of variables that determine the sign (denominators with odd exponent)
	std::deque<Variable *> sign_vars;

	// List of all the explicit variables
	std::deque<ExplicitVariable> all_explicit_vars;
	// List of all the lambda variables
	std::deque<LambdaVariable>   all_lambda_vars;

	// List of error difinitions
	std::deque<ErrorDefinition> all_error_defs;

	std::map<std::string, Variable *> name_2_vars;

	bool is_indirect; // True if function is an indirect predicate
	bool is_lambda;   // True if function defines a lambda

public:
	Predicate(bool _append, bool _output_filtered, bool _output_interval,
	          bool _output_exact, bool _output_expansion);

	bool parseLine(std::ifstream &file, int ln);

	void printErrorBounds();

	void produceAllCode(const std::string &func_name);

protected: /* File and parse **************************************************/
	void openHeader(std::ofstream &header);

	void openFile(std::ofstream &file, const std::string &func_name);

	void parseExplicitVar(std::string &line);

	void parseImplicitVar(std::string &line);

	void parseErrorDefinition(std::string &line);

protected: /* Create codes ****************************************************/
	Variable *getVarByName(const std::string &name);

	std::string createLambdaParamProtoList(const std::string &mytype);

	std::string createParameterProtoList(const std::string &mytype,
	                                     bool               separate_explicit);

	std::string createExpansionProtoList();

	std::string createParameterValueList(bool separate_explicit);

	void produceMultiStageCode(const std::string &func_name,
	                           const std::string &filtered_funcname,
	                           const std::string &interval_funcname,
	                           const std::string &exact_funcname,
	                           const std::string &expansion_funcname,
	                           std::ofstream     &file);

	void producePointArrangement(std::vector<std::vector<size_t>> &arrangements,
	                             std::vector<std::string> &arrangements_str);

	void produceSemiStaticFilter(fpnumber epsilon, int degree,
	                             const std::string &threshold_name,
	                             std::ofstream     &file);

	void produceFilteredCode(const std::string &funcname, std::ofstream &file);

	void produceIntervalCode(const std::string &funcname, std::ofstream &file);

	void produceExpansionCode(const std::string &funcname, std::ofstream &file);

	void produceExactCode(const std::string &funcname, std::ofstream &file);

	void multistagePrototype(const std::string &func_name, std::ofstream &file);

	void filteredPrototype(const std::string &funcname, std::ofstream &file);

	void intervalPrototype(const std::string &funcname, std::ofstream &file);

	void expansionPrototype(const std::string &funcname, std::ofstream &file);

	void exactPrototype(const std::string &funcname, std::ofstream &file);
};