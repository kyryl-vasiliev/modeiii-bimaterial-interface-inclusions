#ifndef INCLUSIONINTEGRALS_H
#define INCLUSIONINTEGRALS_H
 
#include "../../math/templates/Calculation.h"

#include "../../mechanics/templates/BodyPiecewiseHomogeneous.h"
#include "../../mechanics/templates/AnisotropyAngle.h"
#include <complex>
#include <math.h>


#include <optional>

enum class HALFSPACE { HI, LOW };

template<typename T>
std::complex<T> Snzj_plus_ISszj_f5(T t_1, const int &T_index, const T &s_1, const int &vertical_inclusion_index,
	const int &horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T> &b_);// forward declaration

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6(T t_1, const int &T_index, const T &s_1, const int &vertical_inclusion_index,
	const int &horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T> &b_);

template <typename T> //no interface inclusions
std::complex<T> Snzj_plus_ISszj_f5_k_frominclusion(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T> //no interface inclusions
std::complex<T> Snzj_plus_ISszj_f6_k_frominclusion(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_Hi(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_Low(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
std::complex<T> Snzj_plus_ISszj_f6_Hi(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
std::complex<T> Snzj_plus_ISszj_f6_Low(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_frominterface(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
std::complex<T> Snzj_plus_ISszj_f6_frominterface(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_k_frominterface(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_k_frominterface(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_from_interfacial_inclusion(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_from_interfacial_inclusion(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_k_from_interfacial_inclusion(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);
template <typename T>
std::complex<T> Snzj_plus_ISszj_f6_k_from_interfacial_inclusion(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
T W_f5_from_interfacial_inclusion(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);

template <typename T>
T W_f6_from_interfacial_inclusion(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_);


template<typename T>
T mainIntegral(int n, T x) {
	if (n == 0) {
		return 0;
	}
	else {
		return Un<T>(n - 1, x);
	}
}

template<typename T>
bool checkIsInInterface(const int& inclusion_index, const BodyPiecewiseHomogeneous<T>& b_, long double err = 0.0001)
{
	if (fabsl(b_.inclusions[inclusion_index].center.y) < err)
		return true;
	return false;
}

template<typename T>
T jumpFunctionF(int k, T s, int incl_indx, const BodyPiecewiseHomogeneous<T>& b_)
{
	T res = T(0);
	const T& a = b_.inclusions[incl_indx].a;
	s = s / a;
	res = a * Tn<T>(k, s) / sqrtl(1 - s * s);  //f5(s) = sum aj tj(s/al)/sqrt(1 - (s/al))
	return res;
}

template<typename T>
T int_ax(int n, T a, T x) { //int(Tn() ,t=-a..x) x= [-a ..a] a=length of inclusion
	if (n == 0) {
		return a * (T(M_PI) - acosl(x / a));
	}
	else {
		return -a * sinl(n*acosl(x / a)) / n;
	}
}

template<typename T>
T mainDiagonalSyz_f5(int t_index, T s_1,  int main, const BodyPiecewiseHomogeneous<T> &b_) {
	const T &h = b_.inclusions[main].h;
	const T &a = b_.inclusions[main].a;

	const T &a44_incl = b_.inclusions[main].material.a44;
	const T &a45_incl = b_.inclusions[main].material.a45;

	if (!checkIsInInterface(main, b_))
	{

		T res = a45_incl / h / a44_incl * int_ax(t_index, a, s_1 * a); // s [-1..1]
		T snzHiLowf5 = 2.0 * (gaussQuadratureEpsilon(b_.inclusions[main].calc_params.eps, b_.inclusions[main].calc_params.N2_f5,
			Snzj_plus_ISszj_f5_frominterface<T>, t_index, s_1, main, main, b_)).real();
		res += snzHiLowf5;
		return res;
	}
	else //interface inclusion stuff
	{
		AnisotropyAngle<T> ConstsHi, ConstsLow;
		ConstsHi.init(b_.areaHi.material, b_.inclusions[main].phi);
		ConstsLow.init(b_.areaLow.material, b_.inclusions[main].phi);

		T res = a45_incl / h / a44_incl * int_ax(t_index, a, s_1 * a); // s [-1..1]
		// no interface for interface inclusion!
		T chs = (ConstsHi.get_alphaj() * ConstsHi.get_a55j() - ConstsLow.get_alphaj() * ConstsLow.get_a55j());
		T zn = (ConstsHi.get_alphaj() * ConstsHi.get_a55j() + ConstsLow.get_alphaj() * ConstsLow.get_a55j());

		T snzplPlsnzmnf5 =  chs/zn * jumpFunctionF<T>(t_index, s_1 * a, main, b_);
		res += snzplPlsnzmnf5;
		return res;
	}

}

template<typename T>
T nonSymmetricalLoadSyz_f5(T s_1,  int main, const BodyPiecewiseHomogeneous<T> &b_) {
	const T &h = b_.inclusions[main].h;
	const T &a = b_.inclusions[main].a;

	const T &a44_incl = b_.inclusions[main].material.a44;
	const T &a45_incl = b_.inclusions[main].material.a45;

	T res = a45_incl / h / a44_incl * (s_1*a + a )*(b_.inclusions[main].tau_minus - b_.inclusions[main].tau_plus ); // int(tau- -tau+ -a ..a*s_1)  s [-1..1]
	return res;
}

template<typename T>
T nonSymmetricalLoadW_f5(T s_1,  int main, const BodyPiecewiseHomogeneous<T> &b_) {
	AnisotropyAngle<T> ConstsHi, ConstsLow;
	const T &h = b_.inclusions[main].h;
	const T &a = b_.inclusions[main].a;
	const T &a44_incl = b_.inclusions[main].material.a44;
	const T &a45_incl = b_.inclusions[main].material.a45;
	const T &a55_incl = b_.inclusions[main].material.a55;

	ConstsHi.init(b_.areaHi.material, b_.inclusions[main].phi);
	ConstsLow.init(b_.areaHi.material, b_.inclusions[main].phi);

	T r0_2 = fabsl(a44_incl * a55_incl - a45_incl * a45_incl);

	T res = - T(1) / h / a44_incl * r0_2 *  (s_1*a + a )*(b_.inclusions[main].tau_minus - b_.inclusions[main].tau_plus );
	return res;
}

template<typename T>
T mainDiagonalSyz_f6(int t_index, T s_1, int main, const BodyPiecewiseHomogeneous<T> &b_) {
	AnisotropyAngle<T> Consts;
	const T &y = b_.inclusions[main].center.y;
	const T &h = b_.inclusions[main].h;
	const T &a = b_.inclusions[main].a;
	const T &a44_incl = b_.inclusions[main].material.a44;

	if (!checkIsInInterface(main, b_))
	{

		if (y > 0)
		{
			Consts.init(b_.areaHi.material, b_.inclusions[main].phi);
		}
		else
		{
			Consts.init(b_.areaLow.material, b_.inclusions[main].phi);
		}
		T res = T(-1.0) / Consts.get_alphaj() / Consts.get_a55j() * mainIntegral(t_index, s_1);
		T snzHiLowf6 = 2.0 * (gaussQuadratureEpsilon(b_.inclusions[main].calc_params.eps, b_.inclusions[main].calc_params.N2_f6,
			Snzj_plus_ISszj_f6_frominterface<T>, t_index, s_1, main, main, b_)).real(); //only from interface because inclusion influence is main integral itself
		res += snzHiLowf6;
		res = res + T(1) / h / a44_incl * int_ax(t_index, a, s_1 * a); // s [-1..1]
		return res;
	}
	else //interface inclusion stuff
	{
		AnisotropyAngle<T> ConstsHi, ConstsLow;
		ConstsHi.init(b_.areaHi.material, b_.inclusions[main].phi);
		ConstsLow.init(b_.areaLow.material, b_.inclusions[main].phi);

		T res = T(-2.0) / (ConstsHi.get_alphaj()*ConstsHi.get_a55j() + ConstsLow.get_alphaj() * ConstsLow.get_a55j()) * mainIntegral(t_index, s_1);
		res = res + T(1) / h / a44_incl * int_ax(t_index, a, s_1 * a); // s [-1..1]
		return res;
	}
}

template<typename T>
T nonMainDiagonalSyz_f5(int t_index, T s, int main, int other, const BodyPiecewiseHomogeneous<T> &b_) {
	if (!checkIsInInterface(main, b_)) // if main is no interface inclusion
	{
		if (!checkIsInInterface(other, b_)) //if other inclusion in - homogeneous space
		{
			T res = 2.0 * (gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f5,
				Snzj_plus_ISszj_f5<T>, t_index, s, main, other, b_)).real(); //includes influence both from inclusion and interface
			return res;
		}
		else //if other - interfacial inclusion
		{
			T res = 2.0 * (gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f5,
				Snzj_plus_ISszj_f5_from_interfacial_inclusion<T>, t_index, s, main, other, b_)).real();
			return res;
		}
	}
	else //if main is interface inclusion
	{
		if (!checkIsInInterface(other, b_)) //if other inclusion in - homogen
		{
			T res =  (gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f5,
				Snzj_plus_ISszj_f5<T>, t_index, s, main, other, b_)).real(); //includes influence both from inclusion and interface
			return res;

		}
		else //if other - interfacial inclusion
		{
			T res = (gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f5,
				Snzj_plus_ISszj_f5_from_interfacial_inclusion<T>, t_index, s, main, other, b_)).real();
			return res;
		}
	}
}

template<typename T>
T nonMainDiagonalSyz_f6(int t_index, T s, int main, int other, const BodyPiecewiseHomogeneous<T> &b_) {
	if (!checkIsInInterface(main, b_)) //main - homogeneous 
	{
		if (!checkIsInInterface(other, b_))
		{
			T res = 2.0 * (gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f6,
				Snzj_plus_ISszj_f6<T>, t_index, s, main, other, b_)).real();
			return res;
		}
		else //if other - interfacial inclusion there is no stress from interface 
		{
			T res = 2.0 * (gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f6,
				Snzj_plus_ISszj_f6_from_interfacial_inclusion<T>, t_index, s, main, other, b_)).real();
			return res;
		}
	}
	else //main - interfacial
	{
		if (!checkIsInInterface(other, b_)) //other- homogeneous
		{
			T res = (gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f6,
				Snzj_plus_ISszj_f6<T>, t_index, s, main, other, b_)).real();
			return res;
		}
		else //other - interfacial
		{
			T res = (gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f6,
				Snzj_plus_ISszj_f6_from_interfacial_inclusion<T>, t_index, s, main, other, b_)).real();
			return res;
		}
	}
}

//second equation

template<typename T>
T mainDiagonalW_f5(int t_index, T s_1, int main, const BodyPiecewiseHomogeneous<T> &b_) {
	AnisotropyAngle<T> Consts;
	const T &y = b_.inclusions[main].center.y;
	const T &h = b_.inclusions[main].h;
	const T &a = b_.inclusions[main].a;
	const T &a44_incl = b_.inclusions[main].material.a44;
	const T &a45_incl = b_.inclusions[main].material.a45;
	const T &a55_incl = b_.inclusions[main].material.a55;
	if (y > 0 )
	{
		Consts.init(b_.areaHi.material, b_.inclusions[main].phi);
	}
	else
	{
		Consts.init(b_.areaLow.material, b_.inclusions[main].phi);
	}
	T r0_2 = fabsl(a44_incl * a55_incl - a45_incl * a45_incl);

	if (!checkIsInInterface(main, b_))
	{

		T res = Consts.get_a55j() * Consts.get_alphaj() * mainIntegral(t_index, s_1);
		res = res - T(1) / h / a44_incl * r0_2 * int_ax(t_index, a, s_1 * a); // s [-1..1]
		//from interface
		std::complex<T> snz_pl_issz_interface = gaussQuadratureEpsilon(b_.inclusions[main].calc_params.eps, b_.inclusions[main].calc_params.N2_f5,
			Snzj_plus_ISszj_f5_frominterface<T>, t_index, s_1, main, main, b_);
		T snz = snz_pl_issz_interface.real();
		T ssz = snz_pl_issz_interface.imag();
		res += 2.0 * (Consts.get_a45j() * snz + Consts.get_a55j() * ssz);
		return res;
	}
	else //if interfacial stuff
	{
		AnisotropyAngle<T> ConstsHi, ConstsLow;
		ConstsHi.init(b_.areaHi.material, b_.inclusions[main].phi);
		ConstsLow.init(b_.areaLow.material, b_.inclusions[main].phi);
		
		T chs = ConstsHi.get_a55j() * ConstsHi.get_alphaj() * ConstsLow.get_a55j() * ConstsLow.get_alphaj();
		T zn = ConstsHi.get_a55j() * ConstsHi.get_alphaj() + ConstsLow.get_a55j() * ConstsLow.get_alphaj();
		T res = 2 * chs / zn * mainIntegral(t_index, s_1);
		res = res - T(1) / h / a44_incl * r0_2 * int_ax(t_index, a, s_1 * a); // s [-1..1]
		return res;
	}
}

template<typename T>
T mainDiagonalW_f6(int t_index, T s_1, int main, const BodyPiecewiseHomogeneous<T> &b_) {
	
	const T &h = b_.inclusions[main].h;
	const T &a = b_.inclusions[main].a;

	const T &a44_incl = b_.inclusions[main].material.a44;
	const T &a45_incl = b_.inclusions[main].material.a45;

	if (!checkIsInInterface(main, b_))
	{
		T res = a45_incl / h / a44_incl * int_ax(t_index, a, s_1 * a); // s [-1..1]

		std::complex<T> snz_pl_issz_interface = gaussQuadratureEpsilon(b_.inclusions[main].calc_params.eps, b_.inclusions[main].calc_params.N2_f6,
			Snzj_plus_ISszj_f6_frominterface<T>, t_index, s_1, main, main, b_);
		AnisotropyAngle<T> Consts;
		const T& y = b_.inclusions[main].center.y;
		if (y > 0)
		{
			Consts.init(b_.areaHi.material, b_.inclusions[main].phi);
		}
		else
		{
			Consts.init(b_.areaLow.material, b_.inclusions[main].phi);
		}
		T snz = snz_pl_issz_interface.real();
		T ssz = snz_pl_issz_interface.imag();
		res += 2.0 * (Consts.get_a45j() * snz + Consts.get_a55j() * ssz);
		return res;
	}
	else //if interfacial stuff
	{
		AnisotropyAngle<T> ConstsHi, ConstsLow;
		ConstsHi.init(b_.areaHi.material, b_.inclusions[main].phi); //phi should be 0
		ConstsLow.init(b_.areaLow.material, b_.inclusions[main].phi); //phi should be 0

		T res = a45_incl / h / a44_incl * int_ax(t_index, a, s_1 * a); // s [-1..1]
		T chs = (ConstsLow.get_alphaj() * ConstsLow.get_a55j() - ConstsHi.get_alphaj() * ConstsHi.get_a55j() );
		T zn = (ConstsHi.get_alphaj() * ConstsHi.get_a55j() + ConstsLow.get_alphaj() * ConstsLow.get_a55j());

		T wplPlwmnf6 = chs / zn * jumpFunctionF<T>(t_index, s_1 * a, main, b_);
		res += wplPlwmnf6;
		return res;
	}
}

template<typename T>
T nonMainDiagonalW_f5(int t_index, T s, int main, int other, const BodyPiecewiseHomogeneous<T> &b_) {
	if (!checkIsInInterface(main, b_)) // main - regular inclusion
	{
		if (!checkIsInInterface(other, b_)) //other -also regular
		{
			AnisotropyAngle<T> Consts;
			const T& y = b_.inclusions[main].center.y;
			if (y > 0)
			{
				Consts.init(b_.areaHi.material, b_.inclusions[main].phi);
			}
			else
			{
				Consts.init(b_.areaLow.material, b_.inclusions[main].phi);
			}
			std::complex<T> snz_pl_issz = gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f5,
				Snzj_plus_ISszj_f5<T>, t_index, s, main, other, b_);
			T snz = snz_pl_issz.real();
			T ssz = snz_pl_issz.imag();
			T res = 2.0 * (Consts.get_a45j() * snz + Consts.get_a55j() * ssz);
			return res;
		}
		else //other - interface incl
		{
			T wf5 = 2.0 * gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f5,
				W_f5_from_interfacial_inclusion<T>, t_index, s, main, other, b_);
			return wf5;
		}
	}
	else // main - interface incl other - regular incl 
		 //different stuff for homogeneous and homogeneous pertubed solutions
	{
		if (!checkIsInInterface(other, b_)) //other - regular
		{
			AnisotropyAngle<T> ConstsHi, ConstsLow;
			const T& y = b_.inclusions[main].center.y;
			ConstsHi.init(b_.areaHi.material, b_.inclusions[main].phi);
			ConstsLow.init(b_.areaLow.material, b_.inclusions[main].phi);
			
			std::complex<T> snz_pl_isszHi = gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f5,
				Snzj_plus_ISszj_f5_Hi<T>, t_index, s, main, other, b_);
			std::complex<T> snz_pl_isszLow = gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f5,
				Snzj_plus_ISszj_f5_Low<T>, t_index, s, main, other, b_);

			T snz = snz_pl_isszHi.real();
			T ssz = snz_pl_isszHi.imag();
			T res = (ConstsHi.get_a45j() * snz + ConstsHi.get_a55j() * ssz);

			snz = snz_pl_isszLow.real();
			ssz = snz_pl_isszLow.imag();
			res += (ConstsLow.get_a45j() * snz + ConstsLow.get_a55j() * ssz);
			return res;
		}
		else //other - also interface
		{
			T wf5 = gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f5,
				W_f5_from_interfacial_inclusion<T>, t_index, s, main, other, b_);
			return wf5;
		}
	}
}

template<typename T>
T nonMainDiagonalW_f6(int t_index, T s, int main, int other, const BodyPiecewiseHomogeneous<T> &b_) {
	if (!checkIsInInterface(main, b_)) // main - regular inclusion
	{
		if (!checkIsInInterface(other, b_)) //other -also regular
		{
			AnisotropyAngle<T> Consts;
			const T& y = b_.inclusions[main].center.y;
			if (y > 0)
			{
				Consts.init(b_.areaHi.material, b_.inclusions[main].phi);
			}
			else
			{
				Consts.init(b_.areaLow.material, b_.inclusions[main].phi);
			}

			std::complex<T> snz_pl_issz = gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f6,
				Snzj_plus_ISszj_f6<T>, t_index, s, main, other, b_);
			T snz = snz_pl_issz.real();
			T ssz = snz_pl_issz.imag();
			T res = 2.0 * (Consts.get_a45j() * snz + Consts.get_a55j() * ssz);
			return res;
		}
		else
		{
			T wf6 = 2.0 * gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f6,
				W_f6_from_interfacial_inclusion<T>, t_index, s, main, other, b_);
			return wf6;
		}
	}
	else // main = interface inclusion
	{
		if (!checkIsInInterface(other, b_)) //other -also regular
		{
			AnisotropyAngle<T> ConstsHi, ConstsLow;
			const T& y = b_.inclusions[main].center.y;
			ConstsHi.init(b_.areaHi.material, b_.inclusions[main].phi);
			ConstsLow.init(b_.areaLow.material, b_.inclusions[main].phi);

			std::complex<T> snz_pl_isszHi = gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f6,
				Snzj_plus_ISszj_f6_Hi<T>, t_index, s, main, other, b_);
			std::complex<T> snz_pl_isszLow = gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f6,
				Snzj_plus_ISszj_f6_Low<T>, t_index, s, main, other, b_);

			T snz = snz_pl_isszHi.real();
			T ssz = snz_pl_isszHi.imag();
			T res = (ConstsHi.get_a45j() * snz + ConstsHi.get_a55j() * ssz);

			snz = snz_pl_isszLow.real();
			ssz = snz_pl_isszLow.imag();
			res += (ConstsLow.get_a45j() * snz + ConstsLow.get_a55j() * ssz);
			return res;
		}
		else //interface interface
		{
			T wf6 = gaussQuadratureEpsilon(b_.inclusions[other].calc_params.eps, b_.inclusions[other].calc_params.N2_f6,
				W_f6_from_interfacial_inclusion<T>, t_index, s, main, other, b_);
			return wf6;
		}
	}
}

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_k_from_interfacial_inclusion(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {

	using namespace std::complex_literals;
	const int& l = vertical_inclusion_index;
	const int& j = horisontal_inclusion_index;
	const int& k = T_index;

	const T& al = b_.inclusions[l].a;
	const T& phil = b_.inclusions[l].phi;
	const T& x0l = b_.inclusions[l].center.x;
	const T& y0l = b_.inclusions[l].center.y;

	const T& aj = b_.inclusions[j].a;
	const T& phij = b_.inclusions[j].phi;
	const T& x0j = b_.inclusions[j].center.x;
	const T& y0j = b_.inclusions[j].center.y;

	const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));

	AnisotropyAngle<T> ConstsHi, ConstsLow, Consts, Consts3_k;
	ConstsHi.init(b_.areaHi.material, phij);
	ConstsLow.init(b_.areaLow.material, phij);
	if (K == HALFSPACE::HI) 
	{
		Consts = ConstsHi;  Consts3_k = ConstsLow;
	}
	else
	{
		Consts = ConstsLow;  Consts3_k = ConstsHi;
	}

	T nl = 0;
	T sl = al * s_1; //s_1 [-1..1]

	std::complex<T> z0l = x0l + I * y0l;
	std::complex<T> z0j = x0j + I * y0j;

	std::complex<T> sj_plus_inj = (sl * std::exp(I * phil) + z0l - z0j) * std::exp(-I * phij);

	T sj = sj_plus_inj.real();
	T nj = sj_plus_inj.imag();

	std::complex<T> zj = sj + (Consts.get_betaj() + I * Consts.get_alphaj()) * nj;

	T t = t_1 * aj; //param integr

	std::complex<T> t5_k = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - zj)); // sum ak*Tk/(t-z)/sqrt(1-t^2)
	std::complex<T> t5_k_conj = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - std::conj(zj)));
	T chys = T(1.0 / 2.0) * Consts3_k.get_a55j() * Consts3_k.get_alphaj();
	T zn = ConstsHi.get_a55j() * ConstsHi.get_alphaj() + ConstsLow.get_a55j() * ConstsLow.get_alphaj();

	std::complex<T> Snzj_plus_Sszj_f5_interf = aj * chys/zn *
		(Consts.get_gpj() * t5_k - Consts.get_gmj() * t5_k_conj);
	return Snzj_plus_Sszj_f5_interf * std::exp(I * (phil - phij));
}

template <typename T>
T W_f5_from_interfacial_inclusion(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {

	if (!checkIsInInterface(vertical_inclusion_index, b_)) // interfacial inclusion for regular one in hi or low halfspace
	{
		const int& l = vertical_inclusion_index;
		const T& y0l = b_.inclusions[l].center.y;
		AnisotropyAngle<T> Consts;
		if (y0l > 0)
		{
			Consts.init(b_.areaHi.material, b_.inclusions[l].phi);
			std::complex<T> snz_pl_issz = Snzj_plus_ISszj_f5_k_from_interfacial_inclusion(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
			T snz = snz_pl_issz.real();
			T ssz = snz_pl_issz.imag();
			T res = (Consts.get_a45j() * snz + Consts.get_a55j() * ssz);
			return res;
		}
		else
		{
			Consts.init(b_.areaLow.material, b_.inclusions[l].phi);
			std::complex<T> snz_pl_issz = Snzj_plus_ISszj_f5_k_from_interfacial_inclusion(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
			T snz = snz_pl_issz.real();
			T ssz = snz_pl_issz.imag();
			T res = (Consts.get_a45j() * snz + Consts.get_a55j() * ssz);
			return res;
		}
	}
	else // interfacial inclusion for interfacial inclusion
	{
		const int& l = vertical_inclusion_index;
		const T& y0l = b_.inclusions[l].center.y;

		AnisotropyAngle<T> ConstsHi, ConstsLow;
		ConstsHi.init(b_.areaHi.material, b_.inclusions[l].phi);
		ConstsLow.init(b_.areaLow.material, b_.inclusions[l].phi);
		std::complex<T> snz_pl_issz_hi = Snzj_plus_ISszj_f5_k_from_interfacial_inclusion(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
			horisontal_inclusion_index, b_);
		std::complex<T> snz_pl_issz_low = Snzj_plus_ISszj_f5_k_from_interfacial_inclusion(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
			horisontal_inclusion_index, b_);
		T snz = snz_pl_issz_hi.real();
		T ssz = snz_pl_issz_hi.imag();
		T res = ConstsHi.get_a45j() * snz + ConstsHi.get_a55j() * ssz;
		snz = snz_pl_issz_low.real();
		ssz = snz_pl_issz_low.imag();
		res += ConstsLow.get_a45j() * snz + ConstsLow.get_a55j() * ssz;
		return res;
	}
}

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_from_interfacial_inclusion(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {
	
	if (!checkIsInInterface(vertical_inclusion_index, b_)) // interfacial inclusion for regular one in hi or low halfspace
	{
		const int& l = vertical_inclusion_index;
		const T& y0l = b_.inclusions[l].center.y;
		if (y0l > 0)
		{
			return Snzj_plus_ISszj_f5_k_from_interfacial_inclusion(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
		}
		else
		{
			return Snzj_plus_ISszj_f5_k_from_interfacial_inclusion(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
		}
	}
	else // interfacial inclusion for interfacial inclusion
	{
		return Snzj_plus_ISszj_f5_k_from_interfacial_inclusion(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
			horisontal_inclusion_index, b_) +
			Snzj_plus_ISszj_f5_k_from_interfacial_inclusion(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
	}
}

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_frominclusion(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {

	using namespace std::complex_literals;
	const int& l = vertical_inclusion_index;
	const int& j = horisontal_inclusion_index;
	const int& k = T_index;

	const T& al = b_.inclusions[l].a;
	const T& phil = b_.inclusions[l].phi;
	const T& x0l = b_.inclusions[l].center.x;
	const T& y0l = b_.inclusions[l].center.y;

	const T& aj = b_.inclusions[j].a;
	const T& phij = b_.inclusions[j].phi;
	const T& x0j = b_.inclusions[j].center.x;
	const T& y0j = b_.inclusions[j].center.y;

	const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));

	AnisotropyAngle<T> Consts;
	if (!checkIsInInterface(vertical_inclusion_index, b_)) // if not interface inclusion
	{
		//if inclusion from one half space make calculation of influence of inclusion else do nothing probably? 
		if ((y0j * y0l) < 0) // if inclusions from different halfspaces 
		{
			return std::complex < T> {0};
		};

		if ((y0j > 0) && (y0l > 0))  //check is it y0j or y0l
		{
			Consts.init(b_.areaHi.material, phij);
		};

		if ((y0j < 0) && (y0l < 0))
		{
			Consts.init(b_.areaLow.material, phij);
		};
	}
	else // main - interface inclusion just return homogeneous solution for inclusion?  because syz1+syz2
	{    // and homogeneous solution is only in one of halfspaces - in other is equal to 0. 
		 //and note that homogenous pertubed solution in both halfspaces is non zero so should do sum!!!
		if (y0j > 0)
		{
			Consts.init(b_.areaHi.material, phij);
		}
		else
		{
			Consts.init(b_.areaLow.material, phij);
		}
	}

	T nl = 0;
	T sl = al * s_1; //s_1 [-1..1]

	std::complex<T> z0l = x0l + I * y0l;
	std::complex<T> z0j = x0j + I * y0j;

	std::complex<T> sj_plus_inj = (sl * std::exp(I * phil) + z0l - z0j) * std::exp(-I * phij);

	T sj = sj_plus_inj.real();
	T nj = sj_plus_inj.imag();

	std::complex<T> zj = sj + (Consts.get_betaj() + I * Consts.get_alphaj()) * nj;

	T t = t_1 * aj; //param integr

	std::complex<T> t5_k = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - zj)); // sum ak*Tk/(t-z)/sqrt(1-t^2)
	std::complex<T> t5_k_conj = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - std::conj(zj)));

	std::complex<T> Snzj_plus_Sszj_f5 = aj * T(1.0 / 4.0) * (Consts.get_gpj() * t5_k - Consts.get_gmj() * t5_k_conj);
	return Snzj_plus_Sszj_f5 * std::exp(I * (phil - phij));
}

template <typename T> //no interface inclusions
std::complex<T> Snzj_plus_ISszj_f5_k_frominclusion(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {

	using namespace std::complex_literals;
	const int& l = vertical_inclusion_index;
	const int& j = horisontal_inclusion_index;
	const int& k = T_index;

	const T& al = b_.inclusions[l].a;
	const T& phil = b_.inclusions[l].phi;
	const T& x0l = b_.inclusions[l].center.x;
	const T& y0l = b_.inclusions[l].center.y;

	const T& aj = b_.inclusions[j].a;
	const T& phij = b_.inclusions[j].phi;
	const T& x0j = b_.inclusions[j].center.x;
	const T& y0j = b_.inclusions[j].center.y;

	const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));

	AnisotropyAngle<T> Consts;
	if (K == HALFSPACE::HI) 
	{
		if (y0j > 0)
		{
			Consts.init(b_.areaHi.material, phij);
		}
		else
		{
			return std::complex < T> {0};//inclusion in different halfspace answer 0
		}
	};
	if (K == HALFSPACE::LOW)
	{
		if (y0j < 0)
		{
			Consts.init(b_.areaLow.material, phij);
		}
		else
		{
			return std::complex < T> {0}; //inclusion in different halfspace answer 0
		}
	};

	T nl = 0;
	T sl = al * s_1; //s_1 [-1..1]

	std::complex<T> z0l = x0l + I * y0l;
	std::complex<T> z0j = x0j + I * y0j;

	std::complex<T> sj_plus_inj = (sl * std::exp(I * phil) + z0l - z0j) * std::exp(-I * phij);

	T sj = sj_plus_inj.real();
	T nj = sj_plus_inj.imag();

	std::complex<T> zj = sj + (Consts.get_betaj() + I * Consts.get_alphaj()) * nj;

	T t = t_1 * aj; //param integr

	std::complex<T> t5_k = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - zj)); // sum ak*Tk/(t-z)/sqrt(1-t^2)
	std::complex<T> t5_k_conj = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - std::conj(zj)));

	std::complex<T> Snzj_plus_Sszj_f5 = aj * T(1.0 / 4.0) * (Consts.get_gpj() * t5_k - Consts.get_gmj() * t5_k_conj);
	return Snzj_plus_Sszj_f5 * std::exp(I * (phil - phij));
}

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_frominterface(T t_1, const int &T_index, const T &s_1, const int &vertical_inclusion_index,
	const int &horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T> &b_) {

	if (!checkIsInInterface(vertical_inclusion_index, b_)) // if main is not interface inclusion
	{
		const int& l = vertical_inclusion_index;
		const T& y0l = b_.inclusions[l].center.y;
		if (y0l > 0) //hi halfspace
		{
			return Snzj_plus_ISszj_f5_k_frominterface(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
		}
		else //formulas for low halfspace y0l < 0
		{
			return Snzj_plus_ISszj_f5_k_frominterface(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
		}
	}
	else //if main - interfacial inclusion. other - regular so syzhi^l + syzlow^l
	{    //because homogeneous is either hi or low halfspaces but pertubed homogeneous both
		return 
			Snzj_plus_ISszj_f5_k_frominterface(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_)+
			Snzj_plus_ISszj_f5_k_frominterface(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
	}
}

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_k_frominterface(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {

	using namespace std::complex_literals;
	const int& l = vertical_inclusion_index;
	const int& j = horisontal_inclusion_index;
	const int& k = T_index;
	const T& al = b_.inclusions[l].a;
	const T& phil = b_.inclusions[l].phi;
	const T& x0l = b_.inclusions[l].center.x;
	const T& y0l = b_.inclusions[l].center.y;

	const T& aj = b_.inclusions[j].a;
	const T& phij = b_.inclusions[j].phi;
	const T& x0j = b_.inclusions[j].center.x;
	const T& y0j = b_.inclusions[j].center.y;

	const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));
	std::complex<T> z0l = x0l + I * y0l;

	T nl = 0;
	T sl = al * s_1; //s_1 [-1..1]

	std::complex<T> x_plus_iy = sl * std::exp(I * phil) + z0l;
	auto x = x_plus_iy.real();
	auto y = x_plus_iy.imag();

	// material consts
	const T& a45_1 = b_.areaHi.material.a45;
	const T& a55_1 = b_.areaHi.material.a55;
	const T& a45_2 = b_.areaLow.material.a45;
	const T& a55_2 = b_.areaLow.material.a55;
	// alpha material init
	AnisotropyAngle<T> ConstsHi, ConstsLow;
	ConstsHi.init(b_.areaHi.material, 0);
	ConstsLow.init(b_.areaLow.material, 0);
	//consts alpha beta from material hi 1 or low 2 no rotation 
	const auto alpha_1 = ConstsHi.get_alphaj();
	const auto beta_1 = ConstsHi.get_betaj();
	const auto alpha_2 = ConstsLow.get_alphaj();
	const auto beta_2 = ConstsLow.get_betaj();
	//inclusion angles in hi 1 and low 2 same because we use horisontal_inclusion_index
	//works only one of them (same)
	const auto phijk_1 = phij;
	const auto phijk_2 = phij;
	//consts from material hi 1 or low 2 rotation on angle
	ConstsHi.init(b_.areaHi.material, phijk_1); //same angle different materials
	ConstsLow.init(b_.areaLow.material, phijk_2);
	const auto a44jk_1 = ConstsHi.get_a44j();
	const auto a55jk_1 = ConstsHi.get_a55j();
	const auto a44jk_2 = ConstsLow.get_a44j();
	const auto a55jk_2 = ConstsLow.get_a55j();
	const auto alphajk_1 = ConstsHi.get_alphaj();
	const auto betajk_1 = ConstsHi.get_betaj();
	const auto alphajk_2 = ConstsLow.get_alphaj();
	const auto betajk_2 = ConstsLow.get_betaj();

	//center of inclusion (same)
	const auto x0jk_1 = x0j;
	const auto y0jk_1 = y0j;
	const auto x0jk_2 = x0j;
	const auto y0jk_2 = y0j;

	std::complex<T>  t0(0, 0);

	T t = t_1 * aj; //param integr
	//main stuff
	//sincos calc - same due to same angles
	const auto sin_phijk_1 = sinl(phijk_1);
	const auto cos_phijk_1 = cosl(phijk_1);
	const auto sin_2phijk_1 = sinl(2 * phijk_1);
	const auto cos_2phijk_1 = cosl(2 * phijk_1);
	const auto sin_phijk_2 = sinl(phijk_2);
	const auto cos_phijk_2 = cosl(phijk_2);
	const auto sin_2phijk_2 = sinl(2 * phijk_2);
	const auto cos_2phijk_2 = cosl(2 * phijk_2);

	T C1jk_1 = (-1.0 + a44jk_1 / a55jk_1) * sin_phijk_1 * cos_phijk_1 - betajk_1 * cos_2phijk_1;
	T C1jk_2 = (-1.0 + a44jk_2 / a55jk_2) * sin_phijk_2 * cos_phijk_2 - betajk_2 * cos_2phijk_2;

	T C2jk_1 = (-1.0 + a44jk_1 / a55jk_1) * sin_phijk_1 * sin_phijk_1 - betajk_1 * sin_2phijk_1 + 1.0;
	T C2jk_2 = (-1.0 + a44jk_2 / a55jk_2) * sin_phijk_2 * sin_phijk_2 - betajk_2 * sin_2phijk_2 + 1.0;

	std::complex<T> G1jk_1 = a45_1 * C2jk_1 + a55_1 * C1jk_1 + std::complex < T>(1.0i) * (-a55_1 * alphajk_1 + C2jk_1 * a55_2 * alpha_2);
	std::complex<T> G1jk_2 = a45_2 * C2jk_2 + a55_2 * C1jk_2 + std::complex < T>(1.0i) * (-a55_2 * alphajk_2 + C2jk_2 * a55_1 * alpha_1);

	std::complex<T> G2jk_1 = a45_2 * C2jk_2 + a55_2 * C1jk_2 + std::complex < T>(1.0i) * (a55_2 * alphajk_2 + C2jk_2 * a55_2 * alpha_2);
	std::complex<T> G2jk_2 = a45_1 * C2jk_1 + a55_1 * C1jk_1 + std::complex < T>(1.0i) * (a55_1 * alphajk_1 + C2jk_1 * a55_1 * alpha_1);

	T T1jk_1 = alphajk_1 * (t * sin_phijk_1 + y0jk_1);
	T T1jk_2 = alphajk_2 * (t * sin_phijk_2 + y0jk_2);

	T T2jk_1 = (-cos_phijk_1 + betajk_1 * sin_phijk_1) * (t + (cos_phijk_1 - betajk_1 * sin_phijk_1) * x0jk_1 +
		(sin_phijk_1 + betajk_1 * cos_phijk_1) * y0jk_1) - alphajk_1 * alphajk_1 * (x0jk_1 * sin_phijk_1 - y0jk_1 * cos_phijk_1) * sin_phijk_1;
	T T2jk_2 = (-cos_phijk_2 + betajk_2 * sin_phijk_2) * (t + (cos_phijk_2 - betajk_2 * sin_phijk_2) * x0jk_2 +
		(sin_phijk_2 + betajk_2 * cos_phijk_2) * y0jk_2) - alphajk_2 * alphajk_2 * (x0jk_2 * sin_phijk_2 - y0jk_2 * cos_phijk_2) * sin_phijk_2;

	T Bjk_1 = (betajk_1 * sin_phijk_1 - cos_phijk_1) * (betajk_1 * sin_phijk_1 - cos_phijk_1) + alphajk_1 * alphajk_1 * sin_phijk_1 * sin_phijk_1;
	T Bjk_2 = (betajk_2 * sin_phijk_2 - cos_phijk_2) * (betajk_2 * sin_phijk_2 - cos_phijk_2) + alphajk_2 * alphajk_2 * sin_phijk_2 * sin_phijk_2;

	std::complex<T> H1jk_1 = y * Bjk_1 * alpha_1 + T1jk_1 + std::complex <T>(1.0i) * (y * beta_1 * Bjk_1 + T2jk_1 + x * Bjk_1);
	std::complex<T> H1jk_2 = y * Bjk_2 * alpha_2 + T1jk_2 + std::complex <T>(1.0i) * (y * beta_2 * Bjk_2 + T2jk_2 + x * Bjk_2);

	std::complex<T> H2jk_1 = y * Bjk_2 * alpha_1 - T1jk_2 + std::complex <T>(1.0i) * (y * beta_1 * Bjk_2 + T2jk_2 + x * Bjk_2);
	std::complex<T> H2jk_2 = y * Bjk_1 * alpha_2 - T1jk_1 + std::complex <T>(1.0i) * (y * beta_2 * Bjk_1 + T2jk_1 + x * Bjk_1);

	std::complex<T> C1 = 1.0 / 4.0 / M_PI / (a55_1 * alpha_1 + a55_2 * alpha_2);
		if (K == HALFSPACE::HI) //hi halfspace
		{
			if (y0j > 0) //influence of hi inclusion
			{
				std::complex<T> hif5f = C1 * (beta_1 + std::complex <T>(1.0i) * (1.0 - alpha_1)) * G1jk_1 / H1jk_1
					- C1 * (beta_1 + std::complex <T>(1.0i) * (1.0 + alpha_1)) * std::conj(G1jk_1) / std::conj(H1jk_1);
				t0 = hif5f;
			}
			else // influence of low inclusion to hi halfspace
			{
				std::complex<T> hif5s = C1 * (beta_1 + std::complex <T>(1.0i) * (1.0 - alpha_1)) * G2jk_1 / H2jk_1
					- C1 * (beta_1 + std::complex <T>(1.0i) * (1.0 + alpha_1)) * std::conj(G2jk_1) / std::conj(H2jk_1);
				t0 = hif5s;
			}
		}
		else //formulas for low halfspace y0l < 0
		{
			if (y0j > 0) //influence of hi inclusion to low halfspace
			{
				std::complex<T> lowf5f = C1 * (beta_2 + std::complex <T>(1.0i) * (1.0 - alpha_2)) * G2jk_2 / H2jk_2
					- C1 * (beta_2 + std::complex <T>(1.0i) * (1.0 + alpha_2)) * std::conj(G2jk_2) / std::conj(H2jk_2);
				t0 = lowf5f;
			}
			else // influence of low inclusion to low halfspace y0j<0
			{
				std::complex<T> lowf5s = C1 * (beta_2 + std::complex <T>(1.0i) * (1.0 - alpha_2)) * G1jk_2 / H1jk_2
					- C1 * (beta_2 + std::complex <T>(1.0i) * (1.0 + alpha_2)) * std::conj(G1jk_2) / std::conj(H1jk_2);
				t0 = lowf5s;
			}
		}
	t0 = t0 * Tn<T>(k, t / aj); // *f5(t/a)
	t0 = aj * t0; //change interval of integral from -a a t0 -1..1 due to d(aj t) = aj dt
	std::complex<T> Snzj_plus_Sszj_f5_interface = t0 * std::exp(phil * I); // rotation of tensor
	return Snzj_plus_Sszj_f5_interface;
}

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {
	return Snzj_plus_ISszj_f5_frominclusion(t_1, T_index, s_1, vertical_inclusion_index, horisontal_inclusion_index, b_) 
		 + Snzj_plus_ISszj_f5_frominterface(t_1, T_index, s_1, vertical_inclusion_index, horisontal_inclusion_index, b_);
}

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_k(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {
	return Snzj_plus_ISszj_f5_k_frominclusion(K, t_1, T_index, s_1, vertical_inclusion_index, horisontal_inclusion_index, b_)
		+ Snzj_plus_ISszj_f5_k_frominterface(K, t_1, T_index, s_1, vertical_inclusion_index, horisontal_inclusion_index, b_);
}

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_Hi(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {
	return Snzj_plus_ISszj_f5_k(HALFSPACE::HI,t_1,T_index,s_1,vertical_inclusion_index,
		horisontal_inclusion_index,b_);
}

template <typename T>
std::complex<T> Snzj_plus_ISszj_f5_Low(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {
	return Snzj_plus_ISszj_f5_k(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
		horisontal_inclusion_index, b_);
}


template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_frominclusion(T t_1, const int &T_index, const T &s_1, const int &vertical_inclusion_index,
	const int &horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T> &b_) {

	using namespace std::complex_literals;
	const int &l = vertical_inclusion_index;
	const int &j = horisontal_inclusion_index;
	const int &k = T_index;

	const T &al = b_.inclusions[l].a;
	const T &phil = b_.inclusions[l].phi;
	const T &x0l = b_.inclusions[l].center.x;
	const T &y0l = b_.inclusions[l].center.y;

	const T &aj = b_.inclusions[j].a;
	const T &phij = b_.inclusions[j].phi;
	const T &x0j = b_.inclusions[j].center.x;
	const T &y0j = b_.inclusions[j].center.y;
	const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));

	AnisotropyAngle<T> Consts;
	//if inclusion from one half space make calculation of influence of inclusion else do nothing probably? 
	if (!checkIsInInterface(vertical_inclusion_index, b_)) // if not interface inclusion
	{

		if ((y0j * y0l) < 0) // if inclusions from different halfspaces 
		{
			return std::complex < T> {0};
		};

		if ((y0j > 0) && (y0l > 0))  //check is it y0j or y0l
		{
			Consts.init(b_.areaHi.material, phij);
		};

		if ((y0j < 0) && (y0l < 0))
		{
			Consts.init(b_.areaLow.material, phij);
		};
	}
	else
	{
		if (y0j < 0)
		{
			Consts.init(b_.areaLow.material, phij);
		}
		else
		{
			Consts.init(b_.areaHi.material, phij);
		}

	}
	T nl = 0;
	T sl = al * s_1; //s_1 [-1..1]

	std::complex<T> z0l = x0l + I*y0l;
	std::complex<T> z0j = x0j + I*y0j;

	std::complex<T> sj_plus_inj = (sl*std::exp(I*phil) + z0l - z0j)*std::exp(-I*phij);

	T sj = sj_plus_inj.real();
	T nj = sj_plus_inj.imag();

	std::complex<T> zj = sj + (Consts.get_betaj() + I *Consts.get_alphaj() )*nj;

	T t = t_1 * aj; //param integr

	std::complex<T> t6_k = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - zj)); // sum ak*Tk/(t-z)/sqrt(1-t^2)
	std::complex<T> t6_k_conj = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - std::conj(zj)));

	std::complex<T> Snzj_plus_Sszj_f6 = aj * I / T(4.0) / Consts.get_a55j() / Consts.get_alphaj()
		*(Consts.get_gpj()*t6_k + Consts.get_gmj()*t6_k_conj);////!!! f6 +
	return Snzj_plus_Sszj_f6 * std::exp(I*(phil - phij));
}

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_k_frominclusion(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {

	using namespace std::complex_literals;
	const int& l = vertical_inclusion_index;
	const int& j = horisontal_inclusion_index;
	const int& k = T_index;

	const T& al = b_.inclusions[l].a;
	const T& phil = b_.inclusions[l].phi;
	const T& x0l = b_.inclusions[l].center.x;
	const T& y0l = b_.inclusions[l].center.y;

	const T& aj = b_.inclusions[j].a;
	const T& phij = b_.inclusions[j].phi;
	const T& x0j = b_.inclusions[j].center.x;
	const T& y0j = b_.inclusions[j].center.y;
	const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));

	AnisotropyAngle<T> Consts;

	if (K == HALFSPACE::HI)
	{
		if (y0j > 0)
		{
			Consts.init(b_.areaHi.material, phij);
		}
		else
		{
			return std::complex<T>{0};
		}
	}

	if (K == HALFSPACE::LOW)
	{
		if (y0j < 0)
		{
			Consts.init(b_.areaLow.material, phij);
		}
		else
		{
			return std::complex<T>{0};
		}
	}

	T nl = 0;
	T sl = al * s_1; //s_1 [-1..1]

	std::complex<T> z0l = x0l + I * y0l;
	std::complex<T> z0j = x0j + I * y0j;

	std::complex<T> sj_plus_inj = (sl * std::exp(I * phil) + z0l - z0j) * std::exp(-I * phij);

	T sj = sj_plus_inj.real();
	T nj = sj_plus_inj.imag();

	std::complex<T> zj = sj + (Consts.get_betaj() + I * Consts.get_alphaj()) * nj;

	T t = t_1 * aj; //param integr

	std::complex<T> t6_k = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - zj)); // sum ak*Tk/(t-z)/sqrt(1-t^2)
	std::complex<T> t6_k_conj = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - std::conj(zj)));

	std::complex<T> Snzj_plus_Sszj_f6 = aj * I / T(4.0) / Consts.get_a55j() / Consts.get_alphaj()
		* (Consts.get_gpj() * t6_k + Consts.get_gmj() * t6_k_conj);////!!! f6 +
	return Snzj_plus_Sszj_f6 * std::exp(I * (phil - phij));
}

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_k_from_interfacial_inclusion(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {

	using namespace std::complex_literals;
	const int& l = vertical_inclusion_index;
	const int& j = horisontal_inclusion_index;
	const int& k = T_index;

	const T& al = b_.inclusions[l].a;
	const T& phil = b_.inclusions[l].phi;
	const T& x0l = b_.inclusions[l].center.x;
	const T& y0l = b_.inclusions[l].center.y;

	const T& aj = b_.inclusions[j].a;
	const T& phij = b_.inclusions[j].phi;
	const T& x0j = b_.inclusions[j].center.x;
	const T& y0j = b_.inclusions[j].center.y;
	const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));

	AnisotropyAngle<T> ConstsHi, ConstsLow, Consts;
	ConstsHi.init(b_.areaHi.material, phij);
	ConstsLow.init(b_.areaLow.material, phij);
	if (K == HALFSPACE::HI)
	{
		Consts = ConstsHi;
	}
	else
	{
		Consts = ConstsLow;
	}

	T nl = 0;
	T sl = al * s_1; //s_1 [-1..1]

	std::complex<T> z0l = x0l + I * y0l;
	std::complex<T> z0j = x0j + I * y0j;

	std::complex<T> sj_plus_inj = (sl * std::exp(I * phil) + z0l - z0j) * std::exp(-I * phij);

	T sj = sj_plus_inj.real();
	T nj = sj_plus_inj.imag();

	std::complex<T> zj = sj + (Consts.get_betaj() + I * Consts.get_alphaj()) * nj;

	T t = t_1 * aj; //param integr

	std::complex<T> t6_k = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - zj)); // sum ak*Tk/(t-z)/sqrt(1-t^2)
	std::complex<T> t6_k_conj = T(M_1_PI) * (Tn<T>(k, t / aj) / (t - std::conj(zj)));

	T zn = ConstsHi.get_a55j() * ConstsHi.get_alphaj() + ConstsLow.get_a55j() * ConstsLow.get_alphaj();
	std::complex<T> Snzj_plus_Sszj_f6_interf = aj * I / T(2.0) / zn
		* (Consts.get_gpj() * t6_k + Consts.get_gmj() * t6_k_conj);////!!! f6 +
	return Snzj_plus_Sszj_f6_interf * std::exp(I * (phil - phij));
}

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_from_interfacial_inclusion(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {
	if (!checkIsInInterface(vertical_inclusion_index, b_)) // interfacial inclusion for regular one in hi or low halfspace
	{
		const int& l = vertical_inclusion_index; //
		const T& y0l = b_.inclusions[l].center.y;
		if (y0l > 0)
		{
			return Snzj_plus_ISszj_f6_k_from_interfacial_inclusion(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
		}
		else
		{
			return Snzj_plus_ISszj_f6_k_from_interfacial_inclusion(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
		}
	}
	else // interfacial inclusion for interfacial inclusion
	{
		return Snzj_plus_ISszj_f6_k_from_interfacial_inclusion(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
			horisontal_inclusion_index, b_) +
			   Snzj_plus_ISszj_f6_k_from_interfacial_inclusion(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
	}
}

template <typename T>
T W_f6_from_interfacial_inclusion(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {

	if (!checkIsInInterface(vertical_inclusion_index, b_)) // interfacial inclusion for regular one in hi or low halfspace
	{
		const int& l = vertical_inclusion_index; //
		const T& y0l = b_.inclusions[l].center.y;
		AnisotropyAngle<T> Consts;
		if (y0l > 0)
		{
			Consts.init(b_.areaHi.material, b_.inclusions[l].phi);
			std::complex<T> snz_pl_issz = Snzj_plus_ISszj_f6_k_from_interfacial_inclusion(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
			T snz = snz_pl_issz.real();
			T ssz = snz_pl_issz.imag();
			T res = (Consts.get_a45j() * snz + Consts.get_a55j() * ssz); //*2?
			return res;
		}
		else
		{
			Consts.init(b_.areaLow.material, b_.inclusions[l].phi);
			std::complex<T> snz_pl_issz = Snzj_plus_ISszj_f6_k_from_interfacial_inclusion(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
			T snz = snz_pl_issz.real();
			T ssz = snz_pl_issz.imag();
			T res = (Consts.get_a45j() * snz + Consts.get_a55j() * ssz);
			return res;
		}
	}
	else // interfacial inclusion for interfacial inclusion
	{
		const int& l = vertical_inclusion_index; //
		const T& y0l = b_.inclusions[l].center.y;

		AnisotropyAngle<T> ConstsHi, ConstsLow;
		ConstsHi.init(b_.areaHi.material, b_.inclusions[l].phi);
		ConstsLow.init(b_.areaLow.material, b_.inclusions[l].phi);
		std::complex<T> snz_pl_issz_hi = Snzj_plus_ISszj_f6_k_from_interfacial_inclusion(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
			horisontal_inclusion_index, b_);
		std::complex<T> snz_pl_issz_low = Snzj_plus_ISszj_f6_k_from_interfacial_inclusion(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
			horisontal_inclusion_index, b_);
		T snz = snz_pl_issz_hi.real();
		T ssz = snz_pl_issz_hi.imag();
		T res = ConstsHi.get_a45j() * snz + ConstsHi.get_a55j() * ssz;
		snz = snz_pl_issz_low.real();
		ssz = snz_pl_issz_low.imag();
		res += ConstsLow.get_a45j() * snz + ConstsLow.get_a55j() * ssz;
		return res;
	}
}

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_frominterface(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {

	if (!checkIsInInterface(vertical_inclusion_index, b_)) // if not interface inclusion
	{
		const int& l = vertical_inclusion_index;
		const T& y0l = b_.inclusions[l].center.y;

		if (y0l > 0) //hi halfspace
		{
			return Snzj_plus_ISszj_f6_k_frominterface(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
		}
		else //formulas for low halfspace
		{
			return Snzj_plus_ISszj_f6_k_frominterface(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
		}
	}
	else // main - interface inclusion other- regular so hi hi or lowlow
	{
		return Snzj_plus_ISszj_f6_k_frominterface(HALFSPACE::HI, t_1, T_index, s_1, vertical_inclusion_index,
			horisontal_inclusion_index, b_)+
			   Snzj_plus_ISszj_f6_k_frominterface(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
				horisontal_inclusion_index, b_);
	}
}

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_k_frominterface(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {

	using namespace std::complex_literals;
	const int& l = vertical_inclusion_index;
	const int& j = horisontal_inclusion_index;
	const int& k = T_index;
	const T& al = b_.inclusions[l].a;
	const T& phil = b_.inclusions[l].phi;
	const T& x0l = b_.inclusions[l].center.x;
	const T& y0l = b_.inclusions[l].center.y;

	const T& aj = b_.inclusions[j].a;
	const T& phij = b_.inclusions[j].phi;
	const T& x0j = b_.inclusions[j].center.x;
	const T& y0j = b_.inclusions[j].center.y;

	const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));
	std::complex<T> z0l = x0l + I * y0l;

	T nl = 0;
	T sl = al * s_1; //s_1 [-1..1]

	std::complex<T> x_plus_iy = sl * std::exp(I * phil) + z0l;
	auto x = x_plus_iy.real();
	auto y = x_plus_iy.imag();

	// material consts
	const T& a45_1 = b_.areaHi.material.a45;
	const T& a55_1 = b_.areaHi.material.a55;
	const T& a45_2 = b_.areaLow.material.a45;
	const T& a55_2 = b_.areaLow.material.a55;
	// alpha material init
	AnisotropyAngle<T> ConstsHi, ConstsLow;
	ConstsHi.init(b_.areaHi.material, 0);
	ConstsLow.init(b_.areaLow.material, 0);
	//consts alpha beta from material hi 1 or low 2 no rotation 
	const auto alpha_1 = ConstsHi.get_alphaj();
	const auto beta_1 = ConstsHi.get_betaj();
	const auto alpha_2 = ConstsLow.get_alphaj();
	const auto beta_2 = ConstsLow.get_betaj();
	//inclusion angles in hi 1 and low 2 same because we use horisontal_inclusion_index
	//works only one of them (same)
	const auto phijk_1 = phij;
	const auto phijk_2 = phij;
	//consts from material hi 1 or low 2 rotation on angle
	ConstsHi.init(b_.areaHi.material, phijk_1); //same angle different materials
	ConstsLow.init(b_.areaLow.material, phijk_2);
	const auto a44jk_1 = ConstsHi.get_a44j();
	const auto a55jk_1 = ConstsHi.get_a55j();
	const auto a44jk_2 = ConstsLow.get_a44j();
	const auto a55jk_2 = ConstsLow.get_a55j();
	const auto alphajk_1 = ConstsHi.get_alphaj();
	const auto betajk_1 = ConstsHi.get_betaj();
	const auto alphajk_2 = ConstsLow.get_alphaj();
	const auto betajk_2 = ConstsLow.get_betaj();

	//center of inclusion (same)
	const auto x0jk_1 = x0j;
	const auto y0jk_1 = y0j;
	const auto x0jk_2 = x0j;
	const auto y0jk_2 = y0j;

	std::complex<T>  t0(0, 0);

	T t = t_1 * aj; //param integr
	//main stuff
	//sincos calc - same due to same angles
	const auto sin_phijk_1 = sinl(phijk_1);
	const auto cos_phijk_1 = cosl(phijk_1);
	const auto sin_2phijk_1 = sinl(2 * phijk_1);
	const auto cos_2phijk_1 = cosl(2 * phijk_1);
	const auto sin_phijk_2 = sinl(phijk_2);
	const auto cos_phijk_2 = cosl(phijk_2);
	const auto sin_2phijk_2 = sinl(2 * phijk_2);
	const auto cos_2phijk_2 = cosl(2 * phijk_2);

	T C1jk_1 = (-1.0 + a44jk_1 / a55jk_1) * sin_phijk_1 * cos_phijk_1 - betajk_1 * cos_2phijk_1;
	T C1jk_2 = (-1.0 + a44jk_2 / a55jk_2) * sin_phijk_2 * cos_phijk_2 - betajk_2 * cos_2phijk_2;

	T C2jk_1 = (-1.0 + a44jk_1 / a55jk_1) * sin_phijk_1 * sin_phijk_1 - betajk_1 * sin_2phijk_1 + 1.0;
	T C2jk_2 = (-1.0 + a44jk_2 / a55jk_2) * sin_phijk_2 * sin_phijk_2 - betajk_2 * sin_2phijk_2 + 1.0;

	std::complex<T> G1jk_1 = a45_1 * C2jk_1 + a55_1 * C1jk_1 + std::complex < T>(1.0i) * (-a55_1 * alphajk_1 + C2jk_1 * a55_2 * alpha_2);
	std::complex<T> G1jk_2 = a45_2 * C2jk_2 + a55_2 * C1jk_2 + std::complex < T>(1.0i) * (-a55_2 * alphajk_2 + C2jk_2 * a55_1 * alpha_1);

	std::complex<T> G2jk_1 = a45_2 * C2jk_2 + a55_2 * C1jk_2 + std::complex < T>(1.0i) * (a55_2 * alphajk_2 + C2jk_2 * a55_2 * alpha_2);
	std::complex<T> G2jk_2 = a45_1 * C2jk_1 + a55_1 * C1jk_1 + std::complex < T>(1.0i) * (a55_1 * alphajk_1 + C2jk_1 * a55_1 * alpha_1);

	T T1jk_1 = alphajk_1 * (t * sin_phijk_1 + y0jk_1);
	T T1jk_2 = alphajk_2 * (t * sin_phijk_2 + y0jk_2);

	T T2jk_1 = (-cos_phijk_1 + betajk_1 * sin_phijk_1) * (t + (cos_phijk_1 - betajk_1 * sin_phijk_1) * x0jk_1 +
		(sin_phijk_1 + betajk_1 * cos_phijk_1) * y0jk_1) - alphajk_1 * alphajk_1 * (x0jk_1 * sin_phijk_1 - y0jk_1 * cos_phijk_1) * sin_phijk_1;
	T T2jk_2 = (-cos_phijk_2 + betajk_2 * sin_phijk_2) * (t + (cos_phijk_2 - betajk_2 * sin_phijk_2) * x0jk_2 +
		(sin_phijk_2 + betajk_2 * cos_phijk_2) * y0jk_2) - alphajk_2 * alphajk_2 * (x0jk_2 * sin_phijk_2 - y0jk_2 * cos_phijk_2) * sin_phijk_2;

	T Bjk_1 = (betajk_1 * sin_phijk_1 - cos_phijk_1) * (betajk_1 * sin_phijk_1 - cos_phijk_1) + alphajk_1 * alphajk_1 * sin_phijk_1 * sin_phijk_1;
	T Bjk_2 = (betajk_2 * sin_phijk_2 - cos_phijk_2) * (betajk_2 * sin_phijk_2 - cos_phijk_2) + alphajk_2 * alphajk_2 * sin_phijk_2 * sin_phijk_2;

	std::complex<T> H1jk_1 = y * Bjk_1 * alpha_1 + T1jk_1 + std::complex <T>(1.0i) * (y * beta_1 * Bjk_1 + T2jk_1 + x * Bjk_1);
	std::complex<T> H1jk_2 = y * Bjk_2 * alpha_2 + T1jk_2 + std::complex <T>(1.0i) * (y * beta_2 * Bjk_2 + T2jk_2 + x * Bjk_2);

	std::complex<T> H2jk_1 = y * Bjk_2 * alpha_1 - T1jk_2 + std::complex <T>(1.0i) * (y * beta_1 * Bjk_2 + T2jk_2 + x * Bjk_2);
	std::complex<T> H2jk_2 = y * Bjk_1 * alpha_2 - T1jk_1 + std::complex <T>(1.0i) * (y * beta_2 * Bjk_1 + T2jk_1 + x * Bjk_1);

	std::complex<T> C1 = std::complex<T>(1.0i) / std::complex<T>(4.0) / std::complex<T>(M_PI) /
		(a55_1 * alpha_1 + a55_2 * alpha_2);

		if (K == HALFSPACE::HI) //hi halfspace
		{
			if (y0j > 0) //influence of hi inclusion
			{
				std::complex<T> hif6f = C1 / a55jk_1 / alphajk_1 * (beta_1 + std::complex<T>(1.0i) * (1.0 - alpha_1)) * G1jk_1 / H1jk_1
					+ C1 / a55jk_1 / alphajk_1 * (beta_1 + std::complex<T>(1.0i) * (1.0 + alpha_1)) * std::conj(G1jk_1) / std::conj(H1jk_1);
				t0 = hif6f;
			}
			else // influence of low inclusion to hi halfspace
			{
				std::complex<T> hif6s = -C1 / a55jk_2 / alphajk_2 * (beta_1 + std::complex<T>(1.0i) * (1.0 - alpha_1)) * G2jk_1 / H2jk_1
					- C1 / a55jk_2 / alphajk_2 * (beta_1 + std::complex<T>(1.0i) * (1.0 + alpha_1)) * std::conj(G2jk_1) / std::conj(H2jk_1);
				t0 = hif6s;
			}
		}
		else //formulas for low halfspace
		{
			if (y0j > 0) //influence of hi inclusion to low halfspace
			{
				std::complex<T> lowf6f = -C1 / a55jk_1 / alphajk_1 * (beta_2 + std::complex<T>(1.0i) * (1.0 - alpha_2)) * G2jk_2 / H2jk_2
					- C1 / a55jk_1 / alphajk_1 * (beta_2 + std::complex<T>(1.0i) * (1.0 + alpha_2)) * std::conj(G2jk_2) / std::conj(H2jk_2);
				t0 = lowf6f;
			}
			else // influence of low inclusion to low halfspace
			{
				std::complex<T> lowf6s = C1 / a55jk_2 / alphajk_2 * (beta_2 + std::complex<T>(1.0i) * (1.0 - alpha_2)) * G1jk_2 / H1jk_2
					+ C1 / a55jk_2 / alphajk_2 * (beta_2 + std::complex<T>(1.0i) * (1.0 + alpha_2)) * std::conj(G1jk_2) / std::conj(H1jk_2);
				t0 = lowf6s;
			}
		}
	t0 = t0 * Tn<T>(k, t / aj); // *f6(t/a)
	t0 = aj * t0; //change interval of integral from -a a t0 -1..1 due to d(aj t) = aj dt
	std::complex<T> Snzj_plus_Sszj_f6_interface = t0 * std::exp(phil * I); // rotation of tensor
	return Snzj_plus_Sszj_f6_interface;
}


template<typename T>
std::complex<T> Snzj_plus_ISszj_f6(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {
	return Snzj_plus_ISszj_f6_frominclusion(t_1, T_index, s_1, vertical_inclusion_index, horisontal_inclusion_index, b_) + 
		   Snzj_plus_ISszj_f6_frominterface(t_1, T_index, s_1, vertical_inclusion_index, horisontal_inclusion_index, b_);
}

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_k(HALFSPACE K, T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {
	return Snzj_plus_ISszj_f6_k_frominclusion(K, t_1, T_index, s_1, vertical_inclusion_index, horisontal_inclusion_index, b_) +
		   Snzj_plus_ISszj_f6_k_frominterface(K, t_1, T_index, s_1, vertical_inclusion_index, horisontal_inclusion_index, b_);
}

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_Hi(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {
	return Snzj_plus_ISszj_f6_k(HALFSPACE::HI,t_1,T_index,s_1,vertical_inclusion_index,
		horisontal_inclusion_index, b_);
}

template<typename T>
std::complex<T> Snzj_plus_ISszj_f6_Low(T t_1, const int& T_index, const T& s_1, const int& vertical_inclusion_index,
	const int& horisontal_inclusion_index, const BodyPiecewiseHomogeneous<T>& b_) {
	return Snzj_plus_ISszj_f6_k(HALFSPACE::LOW, t_1, T_index, s_1, vertical_inclusion_index,
		horisontal_inclusion_index, b_);
}

template<typename T, typename Iterator>
T frj(T t, Iterator itb, Iterator ite)
{
	int T_index = 0;
	return std::accumulate(itb, ite, T(),
		[&T_index, &t](T sum, T x1) { return sum + Tn<T>(T_index++, t) * x1; }); // returns sum of Aj*Ti no division on sqrt (1-t^2)
}

template<typename T>
std::complex<T> Snzj_plus_ISszj(T t_1, const T& x, const T&y, const int& j, const BodyPiecewiseHomogeneous<T>& b_, const std::vector<T>& A_f5f6) {

	using namespace std::complex_literals;

	const T& aj = b_.inclusions[j].a;
	const T& phij = b_.inclusions[j].phi;
	const T& x0j = b_.inclusions[j].center.x;
	const T& y0j = b_.inclusions[j].center.y;
	const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));

	AnisotropyAngle<T> Consts;
	if (y0j > 0)
	{
		Consts.init(b_.areaHi.material, phij);
	}
	else
	{
		Consts.init(b_.areaLow.material, phij);
	};

	std::complex<T> z0j = x0j + I * y0j;

	std::complex<T> sj_plus_inj = (x + I*y - z0j) * std::exp(- I * phij);

	T sj = sj_plus_inj.real();
	T nj = sj_plus_inj.imag();

	std::complex<T> zj = sj + (Consts.get_betaj() + I * Consts.get_alphaj()) * nj;

	T t = t_1 * aj; //param integr
	
	//iterators !!
	const auto begin_f5j = std::begin(A_f5f6) + std::accumulate(std::begin(b_.inclusions), std::begin(b_.inclusions) + j, 0,
		[](int sum, Inclusion<T> incl) {return sum + incl.calc_params.N1_f5 + 1 + incl.calc_params.N1_f6 + 1; });
	const auto begin_f6j = begin_f5j + b_.inclusions[j].calc_params.N1_f5 + 1;
	const auto end_f5j = begin_f6j;
	const auto end_f6j = begin_f6j + b_.inclusions[j].calc_params.N1_f6 + 1;
	//

	std::complex<T> t6_k = T(M_1_PI) * frj(t / aj, begin_f6j, end_f6j) / (t - zj); // sum ak*Tk/(t-z)/sqrt(1-t^2)
	std::complex<T> t6_k_conj = T(M_1_PI) * frj(t / aj, begin_f6j, end_f6j) / (t - std::conj(zj) );
	std::complex<T> Snzj_plus_Sszj_f6 = aj * I / T(4.0) / Consts.get_a55j() / Consts.get_alphaj()
		* (Consts.get_gpj() * t6_k + Consts.get_gmj() * t6_k_conj);////!!! f6 +

	std::complex<T> t5_k = T(M_1_PI) * frj(t / aj, begin_f5j, end_f5j) / (t - zj); // sum ak*Tk/(t-z)/sqrt(1-t^2)
	std::complex<T> t5_k_conj = T(M_1_PI) * frj(t / aj, begin_f5j, end_f5j) / (t - std::conj(zj) );
	std::complex<T> Snzj_plus_Sszj_f5 = aj * T(1.0 / 4.0) * (Consts.get_gpj() * t5_k - Consts.get_gmj() * t5_k_conj);

	return (Snzj_plus_Sszj_f5 + Snzj_plus_Sszj_f6)* std::exp(- I * phij);
}

template<typename T>
std::complex<T> Sum_Snzj_plus_ISszj(T t_1, const T& x, const T& y, const BodyPiecewiseHomogeneous<T>& b_, const std::vector<T>& A_f5f6) {
	std::complex<T> res{};
	for (std::size_t j = 0; j < b_.inclusions.size(); j++)
		res += Snzj_plus_ISszj(t_1, x, y, j, b_, A_f5f6);
	return res;
}

template<typename T>
std::complex<T> Sum_Snzj_plus_ISszj_noK(T t_1, const T& x, const T& y, const size_t& k, const BodyPiecewiseHomogeneous<T>& b_, const std::vector<T>& A_f5f6) {
	std::complex<T> res{};
	for (std::size_t j = 0; j < b_.inclusions.size(); j++)
	{
		if (j != k)
		{
			res += Snzj_plus_ISszj(t_1, x, y, j, b_, A_f5f6);
		}
	}
	return res;
}
template<typename T>
std::optional<std::complex<T> > isInInclusion(const T& x, const T& y, const Inclusion<T>& in_)
{
	const T& x0 = in_.center.x;
	const T& y0 = in_.center.y;
	const T& phi = in_.phi;
	const T& a = in_.a;
	const T& h = in_.h;
	
	const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));
	std::complex<T> z0 = x0 + I * y0;
	std::complex<T> s_plus_in = (x + I * y  - z0) * std::exp(- I * phi);
	T s = s_plus_in.real();
	T n = s_plus_in.imag();
	long double eps = 1e-7;
	if (((s >= (-a - eps)) && (s <= a + eps)) &&
		((n >= (-h - eps)) && (n <= h + eps)))
	{
		return  s_plus_in;
	}
	else
	{
		return {};
	}
}

template<typename T>
std::optional< std::pair<size_t, std::complex<T>> > findInclusionPointIs(const T& x, const T& y, const BodyPiecewiseHomogeneous<T>& b_)
{
	for (size_t i = 0; i < b_.inclusions.size(); i++)
	{
		auto isInincl = isInInclusion(x, y, b_.inclusions[i]);
		if (isInincl)
			return std::make_pair(i, *isInincl);
	}
	return {};
}

template<typename T, typename Iterator>
T int_fk_U(T t, Iterator itb, Iterator ite)
{
	int T_index = -1;
	return std::accumulate(itb, ite, T(),
		[&T_index, &t](T sum, T x1) { return sum + Un<T>(T_index++, t) * x1; }); // returns sum of Aj*U(j-1,x) no division on sqrt (1-t^2)
}

template<typename T, typename Iterator>
T f_k(T t, Iterator itb, Iterator ite)
{
	int T_index = 0;
	return std::accumulate(itb, ite, T(),
		[&T_index, &t](T sum, T x1) { return sum + Tn<T>(T_index++, t) * x1; })/sqrtl(1 - t*t); // returns sum of Aj*U(j-1,x) no division on sqrt (1-t^2)
}
//stress calculation ???
template<typename T>
std::optional < std::complex<T> >Int_Sum_Snzj_plus_ISszj(const T& x, const T& y, const BodyPiecewiseHomogeneous<T>& b_, const std::vector<T>& A_f5f6) {
	std::complex<T> res = b_.areaHi.syz_isxz.getSnz_plus_iSsz_(0);	
	auto isInInclusion = findInclusionPointIs(x, y, b_);
	if (!isInInclusion)// fre area
	{
		res += gaussQuadratureEpsilon(0.001, 4000, Sum_Snzj_plus_ISszj<T>, x, y, b_, A_f5f6);
	}
	else // inclusionareaprobably should calculaten another way - singularityetc
	{
		const size_t &index = (*isInInclusion).first;
		const auto& sIn = (*isInInclusion).second;
		
		auto s = sIn.real();
		const auto& n = sIn.imag();
		
		const auto& h = b_.inclusions[index].h;
		const auto& a = b_.inclusions[index].a;
		const double eps = 1.e-5;

		res += gaussQuadratureEpsilon(0.001, 4000, Sum_Snzj_plus_ISszj_noK<T>, x, y, index, b_, A_f5f6); // calculate all except singular in xoycoordinates

		if (fabsl(fabsl(s) - a) < eps)
		{
			return {}; //n the edge singularity couldnt calculate
		}
		else
		{
			using namespace std::complex_literals;

			const T& phi = b_.inclusions[index].phi;
			const std::complex<T> I = std::complex<T>(std::complex<T>(1.0i));

			AnisotropyAngle<T> ConstsHi;
			ConstsHi.init(b_.areaHi.material, phi);

			//iterators !!
			const auto begin_f5j = std::begin(A_f5f6) + std::accumulate(std::begin(b_.inclusions), std::begin(b_.inclusions) + index, 0,
				[](int sum, Inclusion<T> incl) {return sum + incl.calc_params.N1_f5 + 1 + incl.calc_params.N1_f6 + 1; });
			const auto begin_f6j = begin_f5j + b_.inclusions[index].calc_params.N1_f5 + 1;
			const auto end_f5j = begin_f6j;
			const auto end_f6j = begin_f6j + b_.inclusions[index].calc_params.N1_f6 + 1;
			//

			auto snz = - 0.5/ ConstsHi.get_alphaj()/ ConstsHi.get_a55j() * int_fk_U( s/a, begin_f6j, end_f6j) ; // sum ak*Uk
			auto dw_ds = 0.5 * ConstsHi.get_alphaj() * ConstsHi.get_a55j() * int_fk_U(s / a, begin_f5j, end_f5j); // sum ak*Uk
			
			
			auto f5 = f_k(s / a, begin_f5j, end_f5j);
			auto f6 = f_k(s / a, begin_f6j, end_f6j);
		
			if (n == 0)
			{
				auto ssz = (dw_ds - snz * ConstsHi.get_a45j()) / ConstsHi.get_a55j();
				res += (snz + I*ssz) * exp(- phi);
				return res;
			}
			if (n > 0)
			{
				snz += -0.5 * f5;
				dw_ds += -0.5 * f6;
				auto ssz = (dw_ds - snz * ConstsHi.get_a45j()) / ConstsHi.get_a55j();
				res += (snz + I *ssz) * exp(- phi);
				return res;
			}
			if (n < 0)
			{
				snz += 0.5 * f5;
				dw_ds += 0.5 * f6;
				auto ssz = (dw_ds - snz * ConstsHi.get_a45j()) / ConstsHi.get_a55j();
				res += (snz + I *ssz) * exp(- phi);
				return res;
			}
		}
	}
	return res;
}



#endif