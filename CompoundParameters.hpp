#ifndef _COMPOUND_PARAMETERS_
#define _COMPOUND_PARAMETERS_


struct Function1Variable {
	virtual type operator()(type T) const = 0;
	virtual ~Function1Variable() { }
};


struct Viscosity : public Function1Variable {
	type A;
	type B;
	type C;

	Viscosity() {}
	Viscosity(type A_): A(A_) {}
	
	type operator()(type T) const {
		return A / T; 
	}

};


struct ConstantViscosity : public Function1Variable {

	type A;

	ConstantViscosity() { };
	ConstantViscosity(type A_) : A(A_) {}

	type operator()(type T) const {
		return A;
	}

};

#endif