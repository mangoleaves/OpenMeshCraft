#include "converter.h"

// floating-point number type in generated code
std::string FT      = "double";
// interval number type (IT) in generated code
std::string IT      = "IT";
// exact number type (ET) in generated code
std::string ET      = "ET";
// boolean type in generated code
std::string Boolean = "bool";
// sign type in generated code
std::string Sign    = "Sign";
// API symbol in generated code
std::string API     = "OpenMeshCraft_API";
// point arrangement used by filter in generated code
std::string PntArr  = "PntArr";
// Exit command in generated code
std::string Exit    = "OMC_EXIT(\"Unsopported points arrangement.\")";

void help_info()
{
	std::cout
	  << "USAGE: converter filename.txt [-a] [-ft] [-it] [-ex] -[et]\n"
	     "If '-a', all the code is appended to 'indirect_predicates.cpp'\n"
	     "If '-ft', output filtered function\n"
	     "If '-it', output interval function\n"
	     "If '-ex', output expansion function\n"
	     "If '-et', output exact function\n";
	exit(1);
}

void error(const char *s, int line)
{
	if (line)
	{
		std::cerr << std::format("Line {}: {}\n", line, s);
		exit(0);
	}
	else
		std::cerr << s << std::endl;
	exit(1);
}

void tokenize(const std::string &line, std::vector<std::string> &tokens,
              char separator)
{
	tokens.clear();
	std::stringstream check(line);
	std::string       intermediate;
	while (std::getline(check, intermediate, separator))
		tokens.push_back(intermediate);
}

void create_heading_comment(std::string &s)
{
	// clang-format off
	s +=
	"/****************************************************************************\n"
	"* Indirect predicates for geometric constructions                           *\n"
	"*                                                                           *\n"
	"* Consiglio Nazionale delle Ricerche                                        *\n"
	"* Istituto di Matematica Applicata e Tecnologie Informatiche                *\n"
	"* Sezione di Genova                                                         *\n"
	"* IMATI-GE / CNR                                                            *\n"
	"*                                                                           *\n"
	"* Authors: Marco Attene                                                     *\n"
	"* Copyright(C) 2019: IMATI-GE / CNR                                         *\n"
	"* All rights reserved.                                                      *\n"
	"*                                                                           *\n"
	"* This program is free software; you can redistribute it and/or modify      *\n"
	"* it under the terms of the GNU Lesser General Public License as published  *\n"
	"* by the Free Software Foundation; either version 3 of the License, or (at  *\n"
	"* your option) any later version.                                           *\n"
	"*                                                                           *\n"
	"* This program is distributed in the hope that it will be useful, but       *\n"
	"* WITHOUT ANY WARRANTY; without even the implied warranty of                *\n"
	"* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  *\n"
	"* General Public License for more details.                                  *\n"
	"*                                                                           *\n"
	"* You should have received a copy of the GNU Lesser General Public License  *\n"
	"* along with this program.  If not, see http://www.gnu.org/licenses/.       *\n"
	"*                                                                           *\n"
	"****************************************************************************/\n"
	"\n/****************************************************************************\n"
	"*                                                                           *\n"
	"* Offset predicates for geometric constructions                             *\n"
	"* Contributors: Guo Jia-Peng                                                *\n"
	"*                                                                           *\n"
	"****************************************************************************/\n"
	"\n/* This code was generated automatically.                                  */\n"
	"/* Do not edit unless you exactly know what you are doing!                 */\n\n";
	// clang-format on
}

fpnumber ulp(fpnumber d)
{
	int exp;
	frexp(d, &exp);
	fpnumber u = ldexp(0.5, exp - 52);
	// Comment this out to ignore INTEL's x387 extended precision
	// u += (u / (1 << 11));
	return u;
}

// If first is true, return empty string, otherwise return comma for separating
// other symbols
std::string sepa_comma(bool &first)
{
	if (first)
	{
		first = false;
		return "";
	}
	else
		return ",";
}

LambdaVariable::LambdaVariable(std::string &n)
  : point_type(n)
{
	if (point_type == "implicitPoint2D")
	{
		dim = 2;
	}
	else if (point_type == "implicitPoint3D")
	{
		dim = 3;
	}
	else
	{
		error("unrecognized implicit point type", -1);
	}
}

std::string LambdaVariable::get_type_string() const
{
	auto pointXD = [this](int x)
	{ return std::format("GenericPoint{}T<{}, {}>", x, IT, ET); };

	if (dim == 2)
		return pointXD(2);
	else if (dim == 3)
		return pointXD(3);
	else
	{
		error("Undefined dimension.", 0);
		return "";
	}
}

std::string LambdaVariable::print_filtered() const
{
	std::stringstream variables;
	bool              first = true;
	for (Variable *v : output_pars)
		variables << sepa_comma(first) << v->name;
	return std::format("{}.getFilteredLambda({}, max_var)", point_name,
	                   variables.str());
}

std::string LambdaVariable::print_interval() const
{
	std::stringstream variables;
	bool              first = true;
	for (Variable *v : output_pars)
		variables << sepa_comma(first) << v->name;
	return std::format("{}.getIntervalLambda({})", point_name, variables.str());
}

std::string LambdaVariable::print_expansion() const
{
	std::stringstream variables;
	bool              first = true;
	for (Variable *v : output_pars)
	{
		if (v->is_lambda || v->is_lambda_d)
			variables << std::format("{0} &{1}, {1}_len", sepa_comma(first), v->name);
		else
			variables << std::format("{0} {1}", sepa_comma(first), v->name);
	}
	return std::format("{}.getExpansionLambda({})", point_name, variables.str());
}

std::string LambdaVariable::print_exact() const
{
	std::stringstream variables;
	bool              first = true;
	for (Variable *v : output_pars)
		variables << sepa_comma(first) << v->name;
	return std::format("{}.getExactLambda({})", point_name, variables.str());
}

ExplicitVariable::ExplicitVariable(std::string &n)
  : point_type(n)
{
	if (point_type == "explicitPoint3D")
		dim = 3;
	else if (point_type == "explicitPoint2D")
		dim = 2;
	else
		error("Unsupported exlicit point type.", 0);
}

std::string ExplicitVariable::get_type_string() const
{
	auto pointXD = [this](size_t x)
	{ return std::format("GenericPoint{}T<{}, {}>", x, IT, ET); };

	return pointXD(dim);
}

std::string ExplicitVariable::get_vars() const
{
	std::stringstream variables;
	bool              first = true;
	for (Variable *v : coords)
		variables << sepa_comma(first) << v->name;
	return variables.str();
}

std::string ExplicitVariable::get_coords() const
{
	std::stringstream variables;
	bool              first = true;
	for (Variable *v : coords)
	{
		std::string var_name = v->name;
		size_t      pos      = var_name.find_last_of("xyz");
		var_name = var_name.substr(0, pos) + "." + var_name.substr(pos, 1) + "()";
		variables << sepa_comma(first) << var_name;
	}
	return variables.str();
}

void ErrorDefinition::parseOneErrorDef(ErrorDef &e, std::string &allvarstuff)
{
	std::vector<std::string> tokens;
	tokenize(allvarstuff, tokens, ';');
	e.error_degree   = atoi(tokens[1].c_str());
	e.size           = atoi(tokens[2].c_str());
	e.fp_error_bound = atof(tokens[3].c_str());
	e.fp_value_bound = atof(tokens[4].c_str());
}

// Plain declaration
Variable::Variable(const std::string &s)
  : name(s)
  , is_part_of_explicit(false)
  , is_part_of_implicit(false)
  , is_lambda(name.starts_with("lambda"))
  , is_lambda_d(name.starts_with("lambda_d"))
  , is_beta(name.starts_with("beta"))
  , is_used(false)
  , actual_length(name + "_len")
{
	op1 = nullptr;
	op2 = nullptr;

	clearError();
}

// Plain explicit declaration
Variable::Variable(std::string &s, ExplicitT)
  : name(s)
  , is_part_of_explicit(true)
  , is_part_of_implicit(false)
  , is_lambda(false)
  , is_lambda_d(false)
  , is_beta(false)
  , is_used(false)
  , actual_length(name + "_len")
{
	op1 = nullptr;
	op2 = nullptr;

	clearError();
}

// Plain lambda declaration
Variable::Variable(std::string &s, ImplicitT)
  : name(s)
  , is_part_of_explicit(false)
  , is_part_of_implicit(true)
  , is_lambda(name.starts_with("l") || name.starts_with("d"))
  , is_lambda_d(name.starts_with("d"))
  , is_beta(name.starts_with("b"))
  , is_used(false)
  , actual_length(name + "_len")
{
	op1 = nullptr;
	op2 = nullptr;

	clearError();
}

// Plain assignment
Variable::Variable(std::string &s, Variable *o1)
  : name(s)
  , is_part_of_explicit(false)
  , is_part_of_implicit(false)
  , is_lambda(name.starts_with("lambda"))
  , is_lambda_d(name.starts_with("lambda_d"))
  , is_beta(name.starts_with("beta"))
  , is_used(false)
  , op1(o1)
  , op2(nullptr)
  , op('=')
{
	if (op1 == nullptr)
		error("operand 1 is null.", -1);
	op1->is_used = true;

	clearError();
}

// Binary operation
Variable::Variable(std::string &s, Variable *o1, Variable *o2, char o)
  : name(s)
  , is_part_of_explicit(false)
  , is_part_of_implicit(false)
  , is_lambda(name.starts_with("lambda"))
  , is_lambda_d(name.starts_with("lambda_d"))
  , is_beta(name.starts_with("beta"))
  , is_used(false)
  , op1(o1)
  , op2(o2)
  , op(o)
{
	if (op1 == nullptr)
		error("operand 1 is null.", -1);
	if (op2 == nullptr)
		error("operand 2 is null.", -1);
	op1->is_used = true;
	op2->is_used = true;
}

void Variable::clearError()
{
	is_a_max = false;

	if (is_part_of_explicit)
	{
		error_evaluated = true;
		size            = 1;
		error_degree    = 1;
		fp_value_bound  = 1;
		fp_error_bound  = 0;
	}
	else if (is_part_of_implicit)
	{
		error_evaluated = true;
		if (is_lambda)
		{
			size           = MAX_STATIC_SIZE;
			error_degree   = 1;
			fp_value_bound = 1;
			fp_error_bound = 1;
		}
		else if (is_beta)
		{
			size           = 1;
			error_degree   = 1;
			fp_value_bound = 1;
			fp_error_bound = 0;
		}
		else
			error("neither lambda nor beta", 0);
		// default values (error_degree, fp_value_bound, fp_error_bound) for
		// lambda/beta is not right, just for some calculation (e.g., find max_var,
		// evaluate size in expansion).
	}
	else if (op1 == nullptr)
	{
		// other input variables
		error_evaluated = true;
		size            = 1;
		error_degree    = 1;
		fp_value_bound  = 1;
		fp_error_bound  = 0;
		if (name == "2")
			fp_value_bound = 2;
	}
	else
	{
		// other middle variables
		error_evaluated = false;
		size            = 1;
		error_degree    = 1;
		fp_value_bound  = 1;
		fp_error_bound  = 0;
	}

	if (op1 != nullptr)
		op1->clearError();
	if (op2 != nullptr)
		op2->clearError();
}

void Variable::setError(int _error_degree, int _size, fpnumber _fp_error_bound,
                        fpnumber _fp_value_bound)
{
	error_degree   = _error_degree;
	size           = _size;
	fp_error_bound = _fp_error_bound;
	fp_value_bound = _fp_value_bound;

	error_evaluated = true;
}

void Variable::propagateError()
{
	if (name == "1" || name == "2")
		goto func_end;
	if (op1 == nullptr && op2 == nullptr)
		goto func_end; // input variable, no operands.

	if (op1 && !op1->error_evaluated)
		op1->propagateError();
	if (op2 && !op2->error_evaluated)
		op2->propagateError();

	if (op == '=')
	{
		// plain assignment, copy evaluated error
		size           = op1->size;
		fp_value_bound = op1->fp_value_bound;
		fp_error_bound = op1->fp_error_bound;
		error_degree   = op1->error_degree;
		if (op1->fp_error_bound == 0)
			op1->is_a_max = true;
	}
	else if (op == '+' || op == '-')
	{
		size = op1->size + op2->size;
		if (op == '-' && op1->fp_error_bound == 0 && op2->fp_error_bound == 0)
		{
			// Translation filter.
			// Commented by Macro Attene:
			// This is the same formula used in FPG. I am not sure it is correct.
			// Why is the value bound 1 and not 2? (1 - (-1)) = 2.
			fp_value_bound = 1;
			fp_error_bound = fp_value_bound * (0.5 * FPN_EPSILON);
			error_degree   = 1;
			is_a_max       = true;
		}
		else
		{
			// Regular sum or subtraction
			fp_value_bound = op1->fp_value_bound + op2->fp_value_bound;
			fpnumber u     = 0.5 * ulp(fp_value_bound);
			fp_value_bound += u;
			fp_error_bound = op1->fp_error_bound + op2->fp_error_bound + u;
			error_degree   = std::max(op1->error_degree, op2->error_degree);
			if (op1->fp_error_bound == 0)
				op1->is_a_max = true;
			if (op2->fp_error_bound == 0)
				op2->is_a_max = true;
		}
	}
	else if (op == '*')
	{
		if (op1->name == std::string("2"))
		{
			size           = op2->size;
			fp_value_bound = op1->fp_value_bound * op2->fp_value_bound;
			fp_error_bound = op2->fp_error_bound;
			error_degree   = op2->error_degree;
		}
		else
		{
			size           = 2 * op1->size * op2->size;
			fp_value_bound = op1->fp_value_bound * op2->fp_value_bound;
			fpnumber u     = 0.5 * ulp(fp_value_bound);
			fp_value_bound += u;
			// The formula herebelow is slightly tighter than FPG's.
			// Original as in FPG is commented below.
			fp_error_bound =
			  op1->fp_error_bound * op2->fp_error_bound +
			  op1->fp_error_bound * (op2->fp_value_bound - op2->fp_error_bound) +
			  op2->fp_error_bound * (op1->fp_value_bound - op1->fp_error_bound) + u;
			// fp_error_bound = op1->fp_error_bound * op2->fp_error_bound +
			//               op1->fp_error_bound * op2->fp_value_bound +
			//               op2->fp_error_bound *op1->fp_value_bound + u;
			error_degree = op1->error_degree + op2->error_degree;
		}

		if (op1->fp_error_bound == 0)
			op1->is_a_max = true;
		if (op2->fp_error_bound == 0)
			op2->is_a_max = true;
	}
	else
		error("Unknown operand", 0);

func_end:
	if (size > 2 * MAX_STATIC_SIZE)
		size = 2 * MAX_STATIC_SIZE;
	if (actual_length.empty())
		actual_length = std::to_string(size);

	error_evaluated = true;
}

Predicate::Predicate(bool _append, bool _output_filtered, bool _output_interval,
                     bool _output_exact, bool _output_expansion)
{
	append           = _append;
	output_filtered  = _output_filtered;
	output_interval  = _output_interval;
	output_exact     = _output_exact;
	output_expansion = _output_expansion;

	is_indirect = false;
	is_lambda   = false;

	// "2" is a special variable representing the constant 2
	all_vars.push_back(Variable(std::string("2")));
	// "1" is a special variable representing the constant 1
	all_vars.push_back(Variable(std::string("1")));
}

Variable *Predicate::getVarByName(const std::string &name)
{
	if (name_2_vars.find(name) == name_2_vars.end())
		return nullptr;
	else
		return name_2_vars.at(name);
}

void Predicate::produceAllCode(const std::string &func_name)
{
	std::string filtered_funcname  = func_name + "_filtered";
	std::string interval_funcname  = func_name + "_interval";
	std::string expansion_funcname = func_name + "_expansion";
	std::string exact_funcname     = func_name + "_exact";

	std::ofstream header;
	openHeader(header);

	std::ofstream file;
	openFile(file, func_name);

	// The implicit point does not have error definitions, ignore filtered func.
	if (output_filtered && !all_lambda_vars.empty() && all_error_defs.empty())
		output_filtered = false;

	// generate prototypes in header
	if (output_filtered)
		filteredPrototype(filtered_funcname, header);
	if (output_interval)
		intervalPrototype(interval_funcname, header);
	if (output_exact)
		exactPrototype(exact_funcname, header);
	if (output_expansion)
		expansionPrototype(expansion_funcname, header);

	if (!is_lambda)
		multistagePrototype(func_name, header);

	// generate implementation in source file
	if (output_filtered)
		produceFilteredCode(filtered_funcname, file);
	if (output_interval)
		produceIntervalCode(interval_funcname, file);
	if (output_exact)
		produceExactCode(exact_funcname, file);
	if (output_expansion)
		produceExpansionCode(expansion_funcname, file);

	if (!is_lambda)
		produceMultiStageCode(func_name, filtered_funcname, interval_funcname,
		                      exact_funcname, expansion_funcname, file);

	header.close();
	file.close();
}

/**
 * @brief parameter prototypes for lambda function,
 * @param mytype number types (double, IT, ET).
 */
std::string Predicate::createLambdaParamProtoList(const std::string &mytype)
{
	bool              first = true;
	std::stringstream variables;
	for (size_t i = 2; i < all_vars.size(); i++)
	{
		const Variable &v = all_vars[i];
		// output lambda or beta var
		if (v.is_lambda || v.is_beta)
			variables << std::format("{}{}& {}", sepa_comma(first), mytype, v.name);

		// ordinary var is just input
		if (v.isInput() && v.is_used && v.name != "1" && v.name != "2")
			variables << std::format("{}{} {}", sepa_comma(first), mytype, v.name);
	}
	return variables.str();
}

/**
 * @param mytype type name (e.g., double, IE, ET)
 * @param separate_explicit separate explicit point to coordinates
 * (e.g., a -> ax, ay, az).
 */
std::string Predicate::createParameterProtoList(const std::string &mytype,
                                                bool separate_explicit)
{
	bool              first = true;
	std::stringstream variables;
	// lambda variables =======================================================
	for (const LambdaVariable &l : all_lambda_vars)
	{
		// Example: , const GenericPoint2/3D & p...
		variables << std::format("{}const {}& {}", sepa_comma(first),
		                         l.get_type_string(), l.point_name);
	}
	// explicit variables =====================================================
	for (const ExplicitVariable &ev : all_explicit_vars)
	{
		if (separate_explicit)
		{
			// Example: , NT px, NT py, NT pz...
			for (const Variable *v : ev.coords)
				variables << std::format("{}{} {}", sepa_comma(first), mytype, v->name);
		}
		else
		{
			// Example: , const GenericPoint2/3D & p...
			variables << std::format("{}const {}& {}", sepa_comma(first),
			                         ev.get_type_string(), ev.point_name);
		}
	}
	// other variables ========================================================
	for (const Variable &v : all_vars)
	{
		if (v.isParameter())
		{
			// Example: , NT t...
			variables << std::format("{}{} {}", sepa_comma(first), mytype, v.name);
		}
	}
	return variables.str();
}

std::string Predicate::createExpansionProtoList()
{
	if (is_lambda)
	{
		std::stringstream variables;
		bool              first = true;
		// input explicit variables
		// example: , double px ...
		for (const Variable &v : all_vars)
		{
			if (v.isInput() && v.name != "1" && v.name != "2" && v.is_used &&
			    !v.is_part_of_implicit)
			{
				variables << std::format("{}double {}", sepa_comma(first), v.name);
			}
		}
		// output lambda and beta variables
		for (const Variable &v : all_vars)
		{
			if (v.is_lambda)
			{
				// example: , double** p, int& p_len ...
				variables << std::format("{0}double** {1}, int& {1}_len",
				                         sepa_comma(first), v.name);
			}
			if (v.is_beta)
			{
				// example: , double& px
				variables << std::format("{}double& {}", sepa_comma(first), v.name);
			}
		}
		return variables.str();
	}
	else
	{
		return createParameterProtoList("double", /*separate_explicit*/ true);
	}
}

std::string Predicate::createParameterValueList(bool separate_explicit)
{
	bool              first = true;
	std::stringstream values;
	for (const LambdaVariable &l : all_lambda_vars)
	{
		values << sepa_comma(first) << l.point_name;
	}
	for (const ExplicitVariable &v : all_explicit_vars)
	{
		if (separate_explicit)
			values << sepa_comma(first) << v.get_vars();
		else
			values << sepa_comma(first) << v.get_coords();
	}
	for (const Variable &v : all_vars)
	{
		if (v.isParameter())
		{
			values << sepa_comma(first) << v.name;
		}
	}
	return values.str();
}

void Predicate::produceMultiStageCode(const std::string &func_name,
                                      const std::string &filtered_funcname,
                                      const std::string &interval_funcname,
                                      const std::string &exact_funcname,
                                      const std::string &expansion_funcname,
                                      std::ofstream     &file)
{
	bool worth_a_semistatic_filter = true;
	// ((*(all_vars.end() - 1)).error_degree <= MAX_WORTH_DEGREE);

	std::string comment       = (worth_a_semistatic_filter) ? ("") : ("//");
	std::string inline_ph     = "";
	std::string template_decl = std::format(
	  "template <typename {}, typename {} {}>", IT, ET,
	  output_filtered ? std::string(", bool WithSSFilter") : std::string(""));
	std::string return_type = Sign;

	// create function that accept explicit points in double

	file << std::format(
	  "{} {}\n{} {}({})\n{{\n", inline_ph, template_decl, return_type, func_name,
	  createParameterProtoList(FT, /*separate_explicit*/ true) +
	    (is_indirect && output_filtered ? std::format(", {} arr", PntArr)
	                                    : std::string("")));

	std::string variables = createParameterValueList(/*separate_explicit*/ true);

	file << std::format("{} ret;\n", return_type);

	if (output_filtered)
	{
		file << "if constexpr (WithSSFilter) {\n";
		file << std::format(
		  "{}ret = {}{}({});\n", comment, filtered_funcname,
		  all_lambda_vars.empty() ? "" : std::format("<{},{}>", IT, ET),
		  variables + (is_indirect ? std::string(", arr") : std::string("")));
		file << comment << "if (is_sign_reliable(ret)) return ret;\n";
		file << "}\n";
	}

	if (output_interval)
	{
		file << std::format("ret = {}{}({});\n", interval_funcname,
		                    all_lambda_vars.empty()
		                      ? std::format("<{}>", IT)
		                      : std::format("<{}, {}>", IT, ET),
		                    variables);
		file << "if (is_sign_reliable(ret)) return ret;\n";
	}

	if (output_exact && !output_expansion)
	{
		file << std::format("return {}{}({});\n", exact_funcname,
		                    all_lambda_vars.empty()
		                      ? std::format("<{}>", ET)
		                      : std::format("<{}, {}>", IT, ET),
		                    variables);
	}

	if (output_expansion)
	{
		file << std::format("return {}{}({});\n", expansion_funcname,
		                    all_lambda_vars.empty()
		                      ? std::string("")
		                      : std::format("<{}, {}>", IT, ET),
		                    variables);
	}

	file << "}\n\n";

	// create function that accept explicit points in point type

	if (!all_explicit_vars.empty())
	{
		file << std::format("{} {}\n{} {}({})\n{{\n", inline_ph, template_decl,
		                    return_type, func_name,
		                    createParameterProtoList(FT,
		                                             /*separate_explicit*/ false) +
		                      (is_indirect && output_filtered
		                         ? std::format(", {} arr", PntArr)
		                         : std::string("")));

		variables =
		  createParameterValueList(/*separate_explicit*/ false) +
		  (is_indirect && output_filtered ? std::string(", arr") : std::string(""));

		file << std::format("return {}{}({});\n", func_name,
		                    std::format("<{}, {} {}>", IT, ET,
		                                output_filtered
		                                  ? std::string(", WithSSFilter")
		                                  : std::string("")),
		                    variables);

		file << "}\n\n";
	}
}

void Predicate::producePointArrangement(
  std::vector<std::vector<size_t>> &arrangements,
  std::vector<std::string>         &arrangements_str)
{
	if (all_lambda_vars.empty() || all_error_defs.empty())
		return;

	size_t lambda_vars_size = all_lambda_vars.size();
	size_t error_defs_size  = all_error_defs.size();
	// assume that error defs are sorted by type_id, e.g., {S(2), L(3), T(4)}

	// arrangement (e.g., SSS, SSL, SST, SLL, SLT, STT, LLL, LLT, LTT, TTT)
	//                    222  223  224  233  234  244  333  334  344  444
	// [highest bit ... lowest bit]
	std::vector<size_t> arr(lambda_vars_size, 0);

	// increase arrangement by one.
	auto inc_arr = [&]() -> bool
	{
		arr[lambda_vars_size - 1] += 1;
		// check if need carry lowest bit.
		if (arr[lambda_vars_size - 1] < error_defs_size)
		{
			return true; // does not need carry, go on
		}
		// need carry the lowest bit, further check if need carry other bits
		size_t i;
		for (i = lambda_vars_size - 1; i > 0 && arr[i] == error_defs_size; --i)
		{
			arr[i] = 0;
			arr[i - 1] += 1;
		}
		// if we need to carry to the highest bit but it can't be carried, stop
		if (i == 0 && arr[0] == error_defs_size)
			return false;
		// adjust lower bits to be not less than current bit
		for (size_t j = i + 1; j < lambda_vars_size; ++j)
			arr[j] = arr[i];
		return true; // go on
	};

	do
	{
		// generate arrangement short names	(e.g., "SSL")
		std::string arr_str = "";
		for (size_t i = 0; i < lambda_vars_size; ++i)
			arr_str += all_error_defs[arr[i]].short_name;
		arrangements.push_back(arr);
		arrangements_str.push_back(arr_str);
	} while (inc_arr());
}

void Predicate::produceSemiStaticFilter(fpnumber epsilon, int degree,
                                        const std::string &threshold_name,
                                        std::ofstream     &file)
{
	int td  = 1;
	int l2d = static_cast<int>(std::floor(std::log2(degree)));
	for (int i = 0; i < l2d; i++)
	{
		file << std::format("{0:} *= {0:};\n", threshold_name);
		td += td;
	}
	for (int i = 0; i < (degree - td); i++)
		file << std::format("{} *= max_var;\n", threshold_name);
	file << std::format("{} *= {};\n", threshold_name, epsilon);
}

void Predicate::produceFilteredCode(const std::string &funcname,
                                    std::ofstream     &file)
{
	// Declare function =======================================================

	std::string inline_ph =
	  (!is_lambda && all_lambda_vars.empty()) ? "inline" : "";
	std::string template_decl =
	  (is_lambda || all_lambda_vars.empty())
	    ? ""
	    : std::format("template <typename {}, typename {}>", IT, ET);
	std::string return_type = is_lambda ? Boolean : Sign;
	std::string variables =
	  is_lambda
	    ? createLambdaParamProtoList(FT) + std::format(", {}& max_var", FT)
	    : createParameterProtoList(FT, /*separate_explicit*/ true) +
	        (is_indirect ? std::format(", {} arr", PntArr) : std::string(""));

	file << std::format("{} {}\n{} {}({})\n{{\n", inline_ph, template_decl,
	                    return_type, funcname, variables);

	// Calculate lambdas ======================================================

	if (is_indirect)
	{
		// declare lambda variables
		// -- declare type
		bool first;
		file << FT;
		// -- declare variables
		first = true;
		for (const Variable &v : all_vars)
			if (v.isInput() && v.is_part_of_implicit)
				file << sepa_comma(first) << v.name;
		file << ", max_var = 0;\n";
		// -- get variables	and check validity
		file << "if (\n";
		first = true;
		for (const LambdaVariable &l : all_lambda_vars)
		{
			file << std::format("{} {}\n", first ? "!" : "|| !", l.print_filtered());
			first = false;
		}
		// -- if invalid, return
		file << ") return Sign::UNCERTAIN;\n\n";
	}

	// Function body - operations =============================================

	for (size_t i = 2; i < all_vars.size(); i++)
	{
		const Variable &v = all_vars[i];
		if (v.isInput())
			continue;

		// -- (type) variable
		if (is_lambda && (v.is_lambda || v.is_beta))
			file << std::format("{}", v.name);
		else
			file << std::format("{} {}", FT, v.name);

		// -- assignment
		if (v.op2 == nullptr)
			file << std::format(" = {};\n", (*v.op1).name);
		else
			file << std::format(" = {} {} {};\n", (*v.op1).name, v.op, (*v.op2).name);
	}

	// Function body - filter calculation =====================================

	// -- find out the output variable
	//    if this is a lambda function, the output var is lambda_d.
	//    otherwise the output var is the last var in assignments.
	std::string outvar_name;
	Variable   *outvar   = nullptr;
	std::string eps_name = "epsilon";
	for (size_t i = 2; i < all_vars.size(); i++)
	{
		outvar      = &all_vars[i];
		outvar_name = all_vars[i].name;
		// if encounter a lambda_d variable...
		// ...and it is not the case "lambda_d = 1"
		if (is_lambda && outvar->is_lambda_d &&
		    (outvar->op2 != nullptr || outvar->op1->name != "1"))
		{
			// we rename the eps_name
			eps_name = outvar_name + "_eps";
			break;
		}
	}

	if (is_lambda && eps_name == "epsilon")
	{
		// this should happen only when this is a lambda function and lambda_d = 1
		file << "\nreturn true;\n}\n\n";
		return;
	}
	// -- clear and propagate to find all max var
	for (Variable &v_ : all_vars)
		v_.clearError();
	// NOTE do not remove propagateError below here, it is used to evaluate other
	// properties of variables except error.
	// for (Variable &v_ : all_vars)
	// 	v_.propagateError();
	outvar->propagateError();

	// -- find out the maximal variable (inputs or difference between inputs)

	bool maxes = std::any_of(all_vars.begin() + 2, all_vars.end(),
	                         [](const Variable &var) { return var.is_a_max; });

	if (is_lambda)
		file << "\ndouble _tmp_fabs;\n";
	else
	{
		if (maxes)
			file << "\ndouble _tmp_fabs;\n";
		if (!is_indirect)
			file << "\ndouble max_var = 0.0;\n";
	}

	for (size_t i = 2; i < all_vars.size(); i++)
		if (all_vars[i].is_a_max)
			file << std::format(
			  "if ((_tmp_fabs = fabs( {} )) > max_var) max_var = _tmp_fabs;\n",
			  all_vars[i].name);

	// -- evaluate error

	if (is_indirect)
	{
		std::vector<std::vector<size_t>> arrangements;
		std::vector<std::string>         arrangements_str;
		producePointArrangement(arrangements, arrangements_str);

		size_t lambda_vars_size = all_lambda_vars.size();

		file << std::format("double {} = max_var;\n", eps_name);
		file << "switch (arr) {\n";
		for (size_t i = 0; i < arrangements.size(); i++)
		{
			const std::vector<size_t> &arr     = arrangements[i];
			const std::string         &arr_str = arrangements_str[i];
			// generate error bound and degree
			// -- clear error
			outvar->clearError();
			// -- set new error to lambda variables
			for (size_t j = 0; j < lambda_vars_size; j++)
			{
				ErrorDefinition &ed = all_error_defs[arr[j]];
				LambdaVariable  &lv = all_lambda_vars[j];
				// -- dim 2: lx,ly,d,bx,by
				// -- dim 3: lx,ly,lz,d,bx,by,bz
				for (size_t k = 0; k < ed.dim * 2 + 1; k++)
				{
					lv.output_pars[k]->setError(ed.pars[k].error_degree, ed.pars[k].size,
					                            ed.pars[k].fp_error_bound,
					                            ed.pars[k].fp_value_bound);
				}
			}
			// -- propagate error
			outvar->propagateError();

			file << std::format("case {}::{}:{{\n", PntArr, arr_str);
			produceSemiStaticFilter(outvar->fp_error_bound, outvar->error_degree,
			                        eps_name, file);
			file << "}\nbreak;\n";

			// TEST
			std::cout << std::format("case {}::{} : FPeps {}  DWeps {}  degree {}\n",
			                         PntArr, arr_str, outvar->fp_error_bound,
			                         outvar->error_degree);
		}
		file << std::format("default:{};}}\n", Exit);
	}
	else // is direct or lambda
	{
		outvar->clearError();
		outvar->propagateError();
		file << std::format("double {} = max_var;\n", eps_name);
		produceSemiStaticFilter(outvar->fp_error_bound, outvar->error_degree,
		                        eps_name, file);
	}

	// Function body - compare with filter and return ==========================

	if (is_lambda)
		file << std::format("\nreturn ({0:} > {0:}_eps || {0:} < -{0:}_eps);\n",
		                    outvar_name);
	else
		file << std::format("\nreturn filter_sign({}, epsilon);\n", outvar_name);

	// Function end
	file << "}\n\n";
}

void Predicate::produceIntervalCode(const std::string &funcname,
                                    std::ofstream     &file)
{
	// Declare function =======================================================

	std::string template_decl =
	  (is_lambda || all_lambda_vars.empty())
	    ? std::format("template <typename {}>", IT)
	    : std::format("template <typename {}, typename {}>", IT, ET);
	std::string inline_ph   = "";
	std::string return_type = is_lambda ? Boolean : Sign;
	std::string variables =
	  is_lambda ? createLambdaParamProtoList(IT)
	            : createParameterProtoList(IT, /*separate_explicit*/ true);

	file << std::format("{} {}\n{} {}({})\n{{\n", inline_ph, template_decl,
	                    return_type, funcname, variables);

	// Calculate lambdas ======================================================

	if (is_indirect)
	{
		// declare lambda variables
		// -- declare type
		file << IT;
		// -- declare variables
		bool first = true;
		for (const Variable &v : all_vars)
			if (v.isInput() && v.is_part_of_implicit)
				file << sepa_comma(first) << v.name;
		file << ";\n";
		// -- get variables	and check validity
		file << "if (\n";
		first = true;
		for (const LambdaVariable &l : all_lambda_vars)
		{
			file << std::format("{} {}\n", first ? "!" : "|| !", l.print_interval());
			first = false;
		}
		// -- if invalid, return
		file << ") return Sign::UNCERTAIN;\n\n";
	}

	// Function body - operations =============================================

	// protect rounding mode in interval calculation
	file << std::format("typename {}::Protector P;\n\n", IT);

	Variable *v = nullptr;
	for (size_t i = 2; i < all_vars.size(); i++)
	{
		v = &all_vars[i];
		if (v->isInput())
			continue;
		std::string  expr;
		std::string &o1 = (*v->op1).name;

		if (v->op == '=' && v->op2 == nullptr)
		{ // assignment
			expr = o1;
		}
		else
		{ // binary operator
			std::string &o2 = (*v->op2).name;

			expr = ((o1 == "2") ? (o2 + " " + v->op + " " + o1)
			                    : (o1 + " " + v->op + " " + o2));
		}
		// check if we need to declare this variable as a new var,
		// then assign to it.
		if (is_lambda && (v->is_lambda || v->is_beta))
			file << std::format("{} = {};\n", v->name, expr);
		else
			file << std::format("{} {} = {};\n", IT, v->name, expr);
	}

	// Function body - dynamic filter =========================================

	if (is_lambda)
	{
		// find "lambda_d" and check if its sign is reliable.
		for (size_t i = 2; i < all_vars.size(); i++)
		{
			const Variable &_v = all_vars[i];
			if (_v.is_lambda_d)
			{
				file << std::format("return {}.is_sign_reliable();\n", _v.name);
				break;
			}
		}
	}
	else
	{
		// check if last variable's sign is reliable
		file << std::format("if (!{}.is_sign_reliable()) return Sign::UNCERTAIN;\n",
		                    v->name);
		// if it is reliable, return its sign.
		file << std::format("return OMC::sign({});\n", v->name);
	}

	// function end ===========================================================
	file << "}\n\n";
}

// TODO: use alloca to allocate memory on stack, avoiding memory on heap.
void Predicate::produceExpansionCode(const std::string &funcname,
                                     std::ofstream     &file)
{
	// Declare function =======================================================

	std::string template_decl =
	  (is_lambda || all_lambda_vars.empty())
	    ? std::string("")
	    : std::format("template <typename {}, typename {}>", IT, ET);
	std::string inline_ph   = is_lambda ? "inline" : "";
	std::string return_type = is_lambda ? "void" : Sign;
	std::string variables   = createExpansionProtoList();

	file << std::format("{} {}\n{} {}({})\n{{\n", inline_ph, template_decl,
	                    return_type, funcname, variables);

	// Calculate stack size ====================================================

	// Account for expansionObject + return_value
	uint32_t fixed_stacksize = 152;
	for (const Variable &v : all_vars)
	{
		if (v.isParameter())
			fixed_stacksize += 8; // Size of a pointer on 64bit systems
	}
	if (is_lambda)
	{
		for (const Variable &v : all_vars)
		{
			if (v.is_lambda)
				fixed_stacksize += 16; // two pointers
			else if (v.is_beta)
				fixed_stacksize += 8; // one reference
		}
	}
	else
	{
		for (const Variable &v : all_vars)
		{
			if (v.isInput() && v.is_part_of_implicit)
			{
				if (v.is_lambda || v.is_lambda_d)
					fixed_stacksize += 16; // two pointers
				else                     //_v.is_beta
					fixed_stacksize += 8;  // one reference
			}
		}
	}

	if (fixed_stacksize >= MAX_STACK_SIZE)
		error("Too many parameters - stack overflow unavoidable.\n", 0);

	uint32_t local_max_size = MAX_STATIC_SIZE * 2;
	uint32_t variable_stacksize;

	// Account for variables

	// -- clear and set default error
	for (Variable &v_ : all_vars)
		v_.clearError();
	for (LambdaVariable &lv : all_lambda_vars)
		for (Variable *v_ : lv.output_pars)
			v_->error_evaluated = true;
	// -- propagate error to get size for each variable
	//    (error degree/bound is not used by expansion code.)
	for (Variable &v_ : all_vars)
		v_.propagateError();

	// -- calculate the local maximal size
	do
	{
		local_max_size /= 2;
		variable_stacksize = 0;
		for (const Variable &v : all_vars)
			if ((!is_lambda || (!v.is_lambda && !v.is_beta)) && v.name != "2" &&
			    v.name != "1")
			{
				if (v.size > static_cast<int>(local_max_size))
					variable_stacksize += (local_max_size + 1) * 8;
				else
					variable_stacksize += (v.size * 8);
			}
			else
				variable_stacksize += 4;
	} while ((fixed_stacksize + variable_stacksize) >= MAX_STACK_SIZE);
	printf("STACK --- %d\n", fixed_stacksize + variable_stacksize);

	// Calculate lambdas ======================================================

	if (is_indirect)
	{
		file << "double return_value = NAN;\n";

#ifdef UNDERFLOW_GUARDING
		file << "#ifdef CHECK_FOR_XYZERFLOWS\n";
		file << "   feclearexcept(FE_ALL_EXCEPT);\n";
		file << "#endif\n";
#endif
		// declare lambda variables and pointers.
		bool first;
		first = true;
		for (const Variable &v : all_vars)
		{
			if (v.isInput() && v.is_part_of_implicit)
			{
				if (v.is_lambda || v.is_lambda_d)
				{
					// if (v.size > local_max_size)
					file << std::format("{0} {1}_p[{2}], *{1} = {1}_p",
					                    (first) ? ("double ") : (", "), v.name,
					                    local_max_size);
					// else
					// file << std::format("{0} {1}_p[{2}], *{1} = {1}_p",
					//                    (first) ? ("double ") : (", "), v.name, v.size);
				}
				else // v.is_beta
				{
					file << std::format("{} {}", (first) ? ("double ") : (", "), v.name);
				}
				first = false;
			}
		}
		file << ";\n";
		// declare lambda variables' length
		first = true;
		for (const Variable &v : all_vars)
		{
			if (v.isInput() && v.is_part_of_implicit)
			{
				if (v.is_lambda || v.is_lambda_d)
				{
					//	if (v.size > local_max_size)
					file << std::format("{} {}_len = {}", (first) ? ("int ") : (", "),
					                    v.name, local_max_size);
					//	else
					//	file << std::format("{} {}_len = {}", (first) ? ("int ") : (", "),
					//                  v.name, v.size);
					first = false;
				}
			}
		}
		file << ";\n";
		// calculate lambda variables.
		for (const LambdaVariable &l : all_lambda_vars)
		{
			file << l.print_expansion() << ";\n";
		}
		// check validity of lambda variables
		file << "if (";
		first = true;
		for (const LambdaVariable &l : all_lambda_vars)
		{
			for (const Variable *_v : l.output_pars)
			{
				if (_v->is_lambda_d)
				{
					file << std::format("{0} ({1}[{1}_len-1] != 0)",
					                    (first) ? ("") : ("&&"), _v->name);
					first = false;
				}
			}
		}
		file << ")\n{\n";
	}

	// seek the first variable that is not input.
	size_t first_op;
	for (first_op = 2; first_op < all_vars.size(); first_op++)
		if (!all_vars[first_op].isInput())
			break;

	// Function body - operations ==============================================

	file << "expansionObject o;\n";

	Variable *v = nullptr;
	size_t    i;

	for (i = 2; i < all_vars.size(); i++)
	{
		v = &all_vars[i];
		if (v->isInput())
			continue;
		// -- get the first operand
		const Variable &op1 = *v->op1;
		std::string o1 = ((is_lambda && op1.is_lambda) ? ("*") : ("")) + op1.name;
		// -- declare variable if it is new.
		int         s1 = op1.size;
		const std::string &al1 = op1.actual_length;
		std::string        lendec;
		if (!is_lambda || (!v->is_lambda && !v->is_beta))
		{
			if (v->size > static_cast<int>(local_max_size))
				file << "double " << v->name << "_p[" << local_max_size << "], *"
				     << v->name << " = " << v->name << "_p;\n";
			else
				file << "double " << v->name << "[" << v->size << "];\n";
			lendec = "int ";
		}
		else
			lendec = "";
		// -- plain assignment
		if (v->op == '=' && v->op2 == nullptr)
		{
			if (o1 == std::string("1"))
			{
				file << "(*" << v->name << ")[0] = 1;\n";
				file << "" << v->name << "_len = 1;\n";
				v->actual_length = "1";
			}
			else if (v->is_beta)
			{
				file << v->name << " = " << v->op1->name << ";\n";
			}
			else
				error("Plain assignment not supported!\n", -1);

			continue;
		}
		// -- binary operator
		const Variable &op2 = *v->op2;
		std::string o2 = ((is_lambda && op2.is_lambda) ? ("*") : ("")) + op2.name;
		int         s2 = op2.size;
		const std::string &al2 = op2.actual_length;

		std::string alb, rlb, lms;
		if (v->is_lambda)
		{
			alb = "";
			rlb = "*";
			lms = v->name + "_len";
		}
		else
		{
			alb = "&";
			rlb = "";
			lms = std::to_string(local_max_size);
		}

		// Known cases
		// [k],n *
		if (o1 == std::string("2"))
		{
			if (v->size > static_cast<int>(local_max_size))
				file << lendec << v->name << "_len = o.Double_With_PreAlloc(" << al2
				     << ", " << o2 << ", " << alb << v->name << ", " << lms << ");\n";
			else
				file << "o.Double(" << al2 << ", " << o2 << ", " << rlb << v->name
				     << ");\n";
			v->actual_length = al2; // v->name + "_len";
		}
		// 1,1 +-*^
		else if (s1 == 1 && s2 == 1)
		{
			if (v->op == '+')
				file << "o.Two_Sum(" << o1 << ", " << o2 << ", " << rlb << v->name
				     << ");\n";
			else if (v->op == '-')
				file << "o.Two_Diff(" << o1 << ", " << o2 << ", " << rlb << v->name
				     << ");\n";
			else if (v->op == '*' && v->op1 != v->op2)
				file << "o.Two_Prod(" << o1 << ", " << o2 << ", " << rlb << v->name
				     << ");\n";
			else if (v->op == '*' && v->op1 == v->op2)
				file << "o.Square(" << o1 << ", " << rlb << v->name << ");\n";
		}
		// 2,1 *
		else if (s1 == 2 && s2 == 1)
		{
			if (v->op == '*')
				file << "o.Two_One_Prod(" << o1 << ", " << o2 << ", " << rlb << v->name
				     << ");\n";
			else if (v->op == '-')
				file << "o.two_One_Diff(" << o1 << ", " << o2 << ", " << rlb << v->name
				     << ");\n";
		}
		// 1,2 *
		else if (s1 == 1 && s2 == 2)
		{
			if (v->op == '*')
				file << "o.Two_One_Prod(" << o2 << ", " << o1 << ", " << rlb << v->name
				     << ");\n";
		}
		// 2,2 +-*^
		else if (s1 == 2 && s2 == 2 && v->op != '*')
		{ // Two_Two_Prod creates length of 8 expansion, which is too long.
			if (v->op == '+')
				file << "o.Two_Two_Sum(" << o1 << ", " << o2 << ", " << rlb << v->name
				     << ");\n";
			else if (v->op == '-')
				file << "o.Two_Two_Diff(" << o1 << ", " << o2 << ", " << rlb << v->name
				     << ");\n";
			else if (v->op == '*')
				file << "o.Two_Two_Prod(" << o1 << ", " << o2 << ", " << rlb << v->name
				     << ");\n";
		}
		// Unknown cases
		// n,n +-*
		else
		{
			if (v->size > static_cast<int>(local_max_size))
			{
				// uint32_t lms = (!is_lambda || !v->is_lambda) ?
				// local_max_size : MAX_STATIC_SIZE;
				if (!is_lambda || !v->is_lambda)
				{
					lms = std::to_string(local_max_size);
				}
				if (v->op == '+')
				{
					if (s2 == 1)
						file << lendec << v->name << "_len = o.Gen_Sum_With_PreAlloc("
						     << al1 << ", " << o1 << ", " << 1 << ", &" << o2 << ", " << alb
						     << v->name << ", " << lms << ");\n";
					else if (s1 == 1)
						file << lendec << v->name << "_len = o.Gen_Sum_With_PreAlloc(" << 1
						     << ", &" << o1 << ", " << al2 << ", " << o2 << ", " << alb
						     << v->name << ", " << lms << ");\n";
					else
						file << lendec << v->name << "_len = o.Gen_Sum_With_PreAlloc("
						     << al1 << ", " << o1 << ", " << al2 << ", " << o2 << ", "
						     << alb << v->name << ", " << lms << ");\n";
				}
				else if (v->op == '-')
				{
					if (s2 == 1)
						file << lendec << v->name << "_len = o.Gen_Diff_With_PreAlloc("
						     << al1 << ", " << o1 << ", " << 1 << ", &" << o2 << ", " << alb
						     << v->name << ", " << lms << ");\n";
					else if (s1 == 1)
						file << lendec << v->name << "_len = o.Gen_Diff_With_PreAlloc(" << 1
						     << ", &" << o1 << ", " << al2 << ", " << o2 << ", " << alb
						     << v->name << ", " << lms << ");\n";
					else
						file << lendec << v->name << "_len = o.Gen_Diff_With_PreAlloc("
						     << al1 << ", " << o1 << ", " << al2 << ", " << o2 << ", "
						     << alb << v->name << ", " << lms << ");\n";
				}
				else if (v->op == '*')
				{
					if (s2 == 1)
						file << lendec << v->name << "_len = o.Gen_Scale_With_PreAlloc("
						     << al1 << ", " << o1 << ", " << o2 << ", " << alb << v->name
						     << ", " << lms << ");\n";
					else if (s1 == 1)
						file << lendec << v->name << "_len = o.Gen_Scale_With_PreAlloc("
						     << al2 << ", " << o2 << ", " << o1 << ", " << alb << v->name
						     << ", " << lms << ");\n";
					else
						file << lendec << v->name << "_len = o.Gen_Product_With_PreAlloc("
						     << al1 << ", " << o1 << ", " << al2 << ", " << o2 << ", "
						     << alb << v->name << ", " << lms << ");\n";
				}
			}
			else if (v->op == '+')
			{
				if (s2 == 1)
					file << lendec << v->name << "_len = o.Gen_Sum(" << al1 << ", " << o1
					     << ", " << 1 << ", &" << o2 << ", " << rlb << v->name << ");\n";
				else if (s1 == 1)
					file << lendec << v->name << "_len = o.Gen_Sum(" << 1 << ", &" << o1
					     << ", " << al2 << ", " << o2 << ", " << rlb << v->name << ");\n";
				else
					file << lendec << v->name << "_len = o.Gen_Sum(" << al1 << ", " << o1
					     << ", " << al2 << ", " << o2 << ", " << rlb << v->name << ");\n";
			}
			else if (v->op == '-')
			{
				if (s2 == 1)
					file << lendec << v->name << "_len = o.Gen_Diff(" << al1 << ", " << o1
					     << ", " << 1 << ", &" << o2 << ", " << rlb << v->name << ");\n";
				else if (s1 == 1)
					file << lendec << v->name << "_len = o.Gen_Diff(" << 1 << ", &" << o1
					     << ", " << al2 << ", " << o2 << ", " << rlb << v->name << ");\n";
				else
					file << lendec << v->name << "_len = o.Gen_Diff(" << al1 << ", " << o1
					     << ", " << al2 << ", " << o2 << ", " << rlb << v->name << ");\n";
			}
			else if (v->op == '*')
			{
				if (s2 == 1)
					file << lendec << v->name << "_len = o.Gen_Scale(" << al1 << ", "
					     << o1 << ", " << o2 << ", " << rlb << v->name << ");\n";
				else if (s1 == 1)
					file << lendec << v->name << "_len = o.Gen_Scale(" << al2 << ", "
					     << o2 << ", " << o1 << ", " << rlb << v->name << ");\n";
				else
					file << lendec << v->name << "_len = o.Gen_Product(" << al1 << ", "
					     << o1 << ", " << al2 << ", " << o2 << ", " << rlb << v->name
					     << ");\n";
			}
			v->actual_length = v->name + "_len";
		}
	}

	file << "\n";

	// Function end ============================================================

	{
		if (!is_lambda)
			file << ((is_indirect) ? ("") : ("double "))
			     << "return_value = " << v->name << "[" << v->actual_length
			     << " - 1];\n";

		// -- free newly allocated memory in function
		for (i--; i >= first_op; i--)
		{
			v = &all_vars[i];
			if (v->size > static_cast<int>(local_max_size) &&
			    (!is_lambda || (!v->is_lambda && !v->is_beta)))
				file << "if (" << v->name << "_p != " << v->name << ") FreeDoubles("
				     << v->name << ");\n";
		}

		if (!is_lambda)
		{
			if (is_indirect)
			{
				file << "}\n\n";
				// -- free allocated memory of lambda variables
				file << std::format("if (!{}::global_cached_values_enabled()){{\n",
				                    all_lambda_vars.front().get_type_string());
				for (const LambdaVariable &lv : all_lambda_vars)
				{
					for (const Variable *_v : lv.output_pars)
					{
						if (_v->is_lambda || _v->is_lambda_d)
						{
							// we will malloc and cache lambda variables,
							// so comment below check and always free.
							// if (_v->size > static_cast<int>(local_max_size))
							file << std::format("if ({0}_p != {0}) FreeDoubles({0});\n",
							                    _v->name);
						}
					}
				}
				file << "}\n\n";
				// -- check result is valid, otherwise call exact function.
#ifdef UNDERFLOW_GUARDING
				file << "#ifdef CHECK_FOR_XYZERFLOWS\n";
				file << "if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW)) "
				        "return ";

				std::string exactname =
				  funcname.substr(0,
				                  funcname.size() - std::string("_expansion").size()) +
				  "_exact";

				std::string template_value = (is_lambda || all_lambda_vars.empty())
				                               ? std::format("<{}>", ET)
				                               : std::format("<{}, {}>", IT, ET);

				file << exactname << template_value << "("
				     << createParameterValueList(/*separate_explicit*/ true) << ");\n";
				file << "#endif\n\n";
#endif
			}

			file << "\nif (return_value > 0) return Sign::POSITIVE;\n";
			file << "if (return_value < 0) return Sign::NEGATIVE;\n";
			file << "if (return_value == 0) return Sign::ZERO;\n";
			// file << "return Sign::UNCERTAIN;\n";
			file << "OMC_EXIT(\"Should not happen.\");\n";
		}
	}

	file << "}\n\n";
}

void Predicate::produceExactCode(const std::string &funcname,
                                 std::ofstream     &file)
{
	// Declare function =====================================================

	std::string template_decl =
	  (is_lambda || all_lambda_vars.empty())
	    ? std::format("template <typename {}>", ET)
	    : std::format("template <typename {}, typename {}>", IT, ET);
	std::string inline_ph   = "";
	std::string return_type = is_lambda ? "void" : Sign;
	std::string variables =
	  is_lambda ? createLambdaParamProtoList(ET)
	            : createParameterProtoList(ET, /*separate_explicit*/ true);

	file << std::format("{} {}\n{} {}({})\n{{\n", inline_ph, template_decl,
	                    return_type, funcname, variables);

	// Calculate lambdas =====================================================

	if (is_indirect)
	{
		// declare lambda variables
		// -- declare type
		file << ET;
		// -- declare variables
		bool first = true;
		for (const Variable &v : all_vars)
			if (v.isInput() && v.is_part_of_implicit)
			{
				file << ((first) ? (" ") : (", ")) << v.name;
				first = false;
			}
		file << ";\n";
		// -- get variables
		for (const LambdaVariable &l : all_lambda_vars)
			file << std::format("{};\n", l.print_exact());
	}

	// Function body - operations =============================================

	Variable *v = nullptr;
	for (size_t i = 2; i < all_vars.size(); i++)
	{
		v = &all_vars[i];
		if (v->isInput())
			continue;
		std::string  expr;
		std::string &o1 = (*v->op1).name;

		if (v->op == '=' && v->op2 == nullptr)
		{ // assignment
			expr = o1;
		}
		else
		{ // binary operator
			std::string &o2 = (*v->op2).name;

			expr = ((o1 == "2") ? (o2 + " " + v->op + " " + o1)
			                    : (o1 + " " + v->op + " " + o2));
		}
		// check if we need to declare this variable as a new var,
		// then assign to it.
		if (is_lambda && (v->is_lambda || v->is_beta))
			file << std::format("{} = {};\n", v->name, expr);
		else
			file << std::format("{} {} = {};\n", ET, v->name, expr);
	}

	// Function body - return	(no filter here) ===============================

	if (!is_lambda)
		file << std::format("return OMC::sign({});\n", v->name);

	// function end ===========================================================
	file << "}\n\n";
}

void Predicate::multistagePrototype(const std::string &func_name,
                                    std::ofstream     &file)
{
	std::string template_decl = std::format(
	  "template <typename {}, typename {} {}>", IT, ET,
	  output_filtered ? std::string(", bool WithSSFilter") : std::string(""));
	std::string return_type = Sign;

	file << std::format(
	  "{} {} {}({});\n\n", template_decl, return_type, func_name,
	  createParameterProtoList(FT, /*separate_explicit*/ true) +
	    (is_indirect && output_filtered ? std::format(", {} arr", PntArr)
	                                    : std::string("")));

	if (!all_explicit_vars.empty())
	{
		file << std::format(
		  "{} {} {}({});\n\n", template_decl, return_type, func_name,
		  createParameterProtoList(FT, /*separate_explicit*/ false) +
		    (is_indirect && output_filtered ? std::format(", {} arr", PntArr)
		                                    : std::string("")));
	}
}

void Predicate::filteredPrototype(const std::string &funcname,
                                  std::ofstream     &file)
{
	// An example:
	// inline template<types...> API return_type func_name(variables...);

	std::string inline_ph =
	  (!is_lambda && all_lambda_vars.empty()) ? "inline" : "";
	std::string template_decl =
	  (is_lambda || all_lambda_vars.empty())
	    ? ""
	    : std::format("template <typename {}, typename {}>", IT, ET);
	std::string export_sign = "";
	std::string return_type = is_lambda ? Boolean : Sign;

	std::string variables =
	  is_lambda
	    ? createLambdaParamProtoList(FT) + std::format(", {}& max_var", FT)
	    : createParameterProtoList(FT, /*separate_explicit*/ true) +
	        (is_indirect ? std::format(", {} arr", PntArr) : std::string(""));

	file << std::format("{} {} {} {} {}({});\n\n", inline_ph, template_decl,
	                    export_sign, return_type, funcname, variables);
}

void Predicate::intervalPrototype(const std::string &funcname,
                                  std::ofstream     &file)
{
	// An example:
	// inline template<type...> API return_type func_name(variables...);

	std::string inline_ph = "";
	std::string template_decl =
	  (is_lambda || all_lambda_vars.empty())
	    ? std::format("template <typename {}>", IT)
	    : std::format("template <typename {}, typename {}>", IT, ET);
	std::string export_sign = "";
	std::string return_type = is_lambda ? Boolean : Sign;

	std::string variables =
	  is_lambda ? createLambdaParamProtoList(IT)
	            : createParameterProtoList(IT, /*separate_explicit*/ true);

	file << std::format("{} {} {} {} {}({});\n\n", inline_ph, template_decl,
	                    export_sign, return_type, funcname, variables);
}

void Predicate::expansionPrototype(const std::string &funcname,
                                   std::ofstream     &file)
{
	// An example:
	// inline template<type...> API return_type func_name(variables...);

	std::string inline_ph =
	  (is_lambda || all_lambda_vars.empty()) ? "inline" : "";
	std::string template_decl =
	  (is_lambda || all_lambda_vars.empty())
	    ? ""
	    : std::format("template <typename {}, typename {}>", IT, ET);
	std::string export_sign = "";
	std::string return_type = is_lambda ? ("void") : Sign;

	std::string variables = createExpansionProtoList();

	file << std::format("{} {} {} {} {}({});\n\n", inline_ph, template_decl,
	                    export_sign, return_type, funcname, variables);
}

void Predicate::exactPrototype(const std::string &funcname, std::ofstream &file)
{
	// An example:
	// inline template<type...> API return_type func_name(variables...);

	std::string inline_ph = "";
	std::string template_decl =
	  (is_lambda || all_lambda_vars.empty())
	    ? std::format("template <typename {}>", ET)
	    : std::format("template <typename {}, typename {}>", IT, ET);
	std::string export_sign = "";
	std::string return_type = is_lambda ? "void" : Sign;

	std::string variables =
	  is_lambda ? createLambdaParamProtoList(ET)
	            : createParameterProtoList(ET, /*separate_explicit*/ true);

	file << std::format("{} {} {} {} {}({});\n\n", inline_ph, template_decl,
	                    export_sign, return_type, funcname, variables);
}

void Predicate::printErrorBounds()
{
	for (size_t i = 2; i < all_vars.size(); i++)
	{
		Variable &v = all_vars[i];
		std::cout << v.name << " " << v.error_degree << " " << v.size << " ";
		std::cout << std::setprecision(std::numeric_limits<fpnumber>::digits10 + 1)
		          << v.fp_error_bound << " ";
		std::cout << std::setprecision(std::numeric_limits<fpnumber>::digits10 + 1)
		          << v.fp_value_bound << "\n";
	}
	std::cout << "NAME ERR_DEGREE EXP_SIZE FP_ERR_BOUND FP_VAL_BOUND\n";
}

void Predicate::openHeader(std::ofstream &header)
{
	std::string heading_comment;
	create_heading_comment(heading_comment);

	bool headerNeedHeadingBlock = true;

	if (append)
	{
		std::ifstream checkopen;
		checkopen.open("indirect_predicates.h");
		if (checkopen.is_open())
		{
			headerNeedHeadingBlock = false;
			checkopen.close();
		}
		header.open("indirect_predicates.h", std::ios_base::app);
	}
	else
	{
		header.open("indirect_predicates.h", std::ios_base::out);
	}

	if (headerNeedHeadingBlock)
	{
		header << heading_comment;
		header << "#include \"implicit_point.h\"\n\n";
	}
}

void Predicate::openFile(std::ofstream &file, const std::string &func_name)
{
	std::string heading_comment;
	create_heading_comment(heading_comment);

	bool fileNeedHeadingBlock = true;
	if (append)
	{
		std::ifstream checkopen;
		checkopen.open("indirect_predicates.hpp");
		if (checkopen.is_open())
		{
			fileNeedHeadingBlock = false;
			checkopen.close();
		}
		file.open("indirect_predicates.hpp", std::ios_base::app);
	}
	else
	{
		file.open(func_name + ".hpp", std::ios_base::out);
	}

	if (fileNeedHeadingBlock)
	{
		file << heading_comment;
		file << "#include \"implicit_point.h\"\n\n";
		file << "#pragma intrinsic(fabs)\n\n";
#ifdef UNDERFLOW_GUARDING
		file
		  << "// Uncomment the following to activate overflow/underflow checks\n";
		file << "#define CHECK_FOR_XYZERFLOWS\n\n";
#endif
	}
}

void Predicate::parseExplicitVar(std::string &line)
{
	// For an explicit var, we expect it is declared in the form:
	// explicitPoint2D(pnt_name:pnt_x,pnt_y) or
	// explicitPoint3D(pnt_name:pnt_x,pnt_y,pnt_z).

	// find "explicitPoint2D" or "explicitPoint3D"
	size_t      pos        = line.find('(');
	std::string pnt_type   = line.substr(0, pos);
	// remaining part goes to "defs"
	std::string defs       = line.substr(pos + 1, line.size() - (pos + 2));
	// find "pnt_name", e.g., p
	pos                    = defs.find(':');
	std::string pnt_name   = defs.substr(0, pos);
	// find "pnt_coords", e.g., px,py,pz
	std::string pnt_coords = defs.substr(pos + 1, defs.size() - (pos + 1));

	// put this new point to the global list
	all_explicit_vars.push_back(ExplicitVariable(pnt_type));
	ExplicitVariable &v = all_explicit_vars.back();

	// get new point's name
	std::vector<std::string> tokens;
	tokenize(pnt_name, tokens, ',');
	v.point_name = tokens[0];

	// store coordinates as variables
	if (pnt_coords.find_first_of(',') != std::string::npos)
	{ // multiple coordinates
		tokenize(pnt_coords, tokens, ',');

		for (size_t i = 0; i < tokens.size(); i++)
		{
			all_vars.push_back(Variable(tokens[i], Variable::ExplicitT()));
			v.coords.push_back(&all_vars.back());
			name_2_vars[all_vars.back().name] = &all_vars.back();
		}
	}
	else
	{ // single coordinate
		all_vars.push_back(Variable(pnt_coords, Variable::ExplicitT()));
		v.coords.push_back(&all_vars.back());
		name_2_vars[all_vars.back().name] = &all_vars.back();
	}
}

void Predicate::parseImplicitVar(std::string &line)
{
	size_t      pos        = line.find('(');
	// implicitPoint2D or implicitPoint3D
	std::string pnt_type   = line.substr(0, pos);
	// lambda variable is seen as a function
	std::string parameters = line.substr(pos + 1, line.size() - (pos + 2));
	pos                    = parameters.find(':');
	// input point name, e.g., p1
	std::string input_pars = parameters.substr(0, pos);
	// output implicit values, e.g., l1x,l1y,l1z,d1,b1x,b1y,b1z
	std::string output_pars =
	  parameters.substr(pos + 1, parameters.size() - (pos + 1));

	all_lambda_vars.push_back(LambdaVariable(pnt_type));
	LambdaVariable &l = all_lambda_vars.back();

	std::vector<std::string> tokens;
	tokenize(input_pars, tokens, ',');

	l.point_name = tokens[0];

	if (output_pars.find_first_of(',') != std::string::npos)
	{
		tokenize(output_pars, tokens, ',');

		for (size_t i = 0; i < tokens.size(); i++)
		{
			all_vars.push_back(Variable(tokens[i], Variable::ImplicitT()));
			l.output_pars.push_back(&all_vars.back());
			name_2_vars[all_vars.back().name] = &all_vars.back();
		}
	}
	else
	{
		error("can't parse lambda variable because fail to tokenize output_pars.",
		      -1);
	}

	is_indirect = true;
}

void Predicate::parseErrorDefinition(std::string &line)
{
	size_t      pos        = line.find('(');
	// lambda variable is seen as a function
	std::string parameters = line.substr(pos + 1, line.size() - (pos + 2));
	pos                    = parameters.find(':');
	// input point name, e.g., SSI, LPI, TPI
	std::string input_pars = parameters.substr(0, pos);
	// output lambda values, e.g., x,y,z,d
	std::string output_pars =
	  parameters.substr(pos + 1, parameters.size() - (pos + 1));

	all_error_defs.emplace_back();
	ErrorDefinition &ed = all_error_defs.back();

	std::vector<std::string> tokens;

	tokenize(input_pars, tokens, ',');

	ed.name       = tokens[0];
	ed.short_name = tokens[1];
	ed.dim        = atoi(tokens[2].c_str());

	tokenize(output_pars, tokens, ';');

	ed.pars.reserve(ed.dim * 2 + 1);

	// parse error definition for lx, ly, lz, d
	for (size_t i = 0; i < tokens.size(); i += 5)
	{
		std::string allvarstuff = tokens[i] + ";" + tokens[i + 1] + ";" +
		                          tokens[i + 2] + ";" + tokens[i + 3] + ";" +
		                          tokens[i + 4];
		ed.pars.emplace_back();
		ed.parseOneErrorDef(ed.pars.back(), allvarstuff);
	}
	// directly add error definition for bx,by,bz
	for (size_t i = 0; i < ed.dim; i++)
	{
		ed.pars.emplace_back();
		ed.pars.back().error_degree   = 1;
		ed.pars.back().size           = 1;
		ed.pars.back().fp_error_bound = 0;
		ed.pars.back().fp_value_bound = 1;
	}

	if (ed.dim == 2)
		PntArr = "PntArr2";
	else if (ed.dim == 3)
		PntArr = "PntArr3";
}

bool Predicate::parseLine(std::ifstream &file, int ln)
{
	std::string              line;
	std::vector<std::string> all_toks, tokens;
	if (!std::getline(file, line))
		return false;
	tokenize(line, all_toks, ' ');

	for (unsigned int i = 0; i < all_toks.size(); i++)
		if (all_toks[i][0] == '/' && all_toks[i][1] == '/')
			break; // meet a comment
		else
			tokens.push_back(all_toks[i]);

	if (tokens.size() == 0)
		return true; // Skip blank lines

	if (tokens.size() == 1)
	{
		size_t pos = tokens[0].find('(');
		if (pos == std::string::npos) // Single parameter
		{
			all_vars.push_back(Variable(tokens[0]));
			name_2_vars[tokens[0]] = &all_vars.back();
		}
		else if (tokens[0].starts_with("explicitPoint"))
		{
			parseExplicitVar(tokens[0]);
		}
		else if (tokens[0].starts_with("implicitPoint"))
		{
			parseImplicitVar(tokens[0]);
		}
		else if (tokens[0].starts_with("errorDefinition"))
		{
			parseErrorDefinition(tokens[0]);
		}
		else
		{
			error("Error: unrecognized var", ln);
		}
	}
	else
	{
		if (tokens[1] == "=") // Assignment
		{
			if (tokens.size() == 5) // Binary op
			{
				std::string &newvarname = tokens[0];
				// CHECK if variable is declared
				if (getVarByName(newvarname) != nullptr)
					error("Error: variable already declared", ln);
				// CHECK if operands is valid
				Variable *op1 = getVarByName(tokens[2]);
				if (op1 == nullptr)
					error("Error: unknown variable used", ln);
				Variable *op2 = getVarByName(tokens[4]);
				if (op2 == nullptr)
					error("Error: unknown variable used", ln);
				// CHECK if operator is valid
				char op = tokens[3][0];
				if (support_operators.find_first_of(op) == std::string::npos)
					error("Error: operator not supported", ln);
				// Put it into the global list
				all_vars.push_back(Variable(newvarname, op1, op2, op));
				name_2_vars[newvarname] = &all_vars.back();
			}
			else if (tokens.size() == 3) // Plain assignment
			{
				std::string &newvarname = tokens[0];
				// CHECK if variable is declared
				if (getVarByName(newvarname) != nullptr)
					error("Error: variable already declared", ln);
				// CHECK if operands is valid
				Variable *op1 = getVarByName(tokens[2]);
				if (op1 == nullptr)
					error("Error: unknown variable used", ln);
				// Put it into the global list
				all_vars.push_back(Variable(newvarname, op1));
				name_2_vars[newvarname] = &all_vars.back();
			}
			else
				error("Error: unexpected number of tokens", ln);
		}
		else // Parameter list
		{
			for (size_t i = 0; i < tokens.size(); i++)
			{
				all_vars.push_back(Variable(tokens[i]));
				name_2_vars[tokens[i]] = &all_vars.back();
			}
		}
	}

	return true;
}