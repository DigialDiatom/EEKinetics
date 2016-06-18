#pragma once
#include <tuple>
namespace ee_kin {
template<int ifactor = 1, int idivisor = 1>
struct Substrate {
	Substrate(float& idestination) : destination(idestination) {}
	float& destination;
	static constexpr float factor = float(ifactor) / float(idivisor);
};

template<int ifactor = 1, int idivisor = 1>
struct Product {
	Product(float& idestination) : destination(idestination) {}
	float& destination;
	static constexpr float factor = float(ifactor) / float(idivisor);
};

template<size_t i, typename T>
constexpr float getFactor() {
  using R = typename std::tuple_element<i, T>::type;
  return R::factor;
}

template <typename Tup, size_t...s>
float multiplyTuple(Tup const &tup, std::index_sequence<s...>)
{
	auto res = 1.0f;
	auto x = { (res *= std::get<s>(tup).destination)... };
	(void)x;
	return res;
}

template <typename Tup>
float multiplyTuple(Tup const &tup)
{
	return multiplyTuple(tup, std::make_index_sequence<std::tuple_size<Tup>::value>());
}

template <typename Tup, size_t...s>
float sumTuple(Tup const &tup, std::index_sequence<s...>)
{
	auto res = 1.0f;
	auto x = { (res += std::get<s>(tup).destination)... };
	(void)x;
	return res;
}

template <typename Tup>
float sumTuple(Tup const &tup)
{
	return sumTuple(tup, std::make_index_sequence<std::tuple_size<Tup>::value>());
}

template <typename Tup, size_t...s>
void decreaseSubstrates(Tup &tup, float rate, std::index_sequence<s...>)
{
	//auto res = 1.0f;
	auto x = { (std::get<s>(tup).destination -= rate * getFactor<s, Tup>())... };
	(void)x;
	//return res;
}

template <typename Tup>
void decreaseSubstrates(Tup &tup, float rate)
{
	decreaseSubstrates(tup, rate, std::make_index_sequence<std::tuple_size<Tup>::value>());
}

template <typename Tup, size_t...s>
void increaseProducts(Tup &tup, float rate, std::index_sequence<s...>)
{
	//auto res = 1.0f;
	auto x = { (std::get<s>(tup).destination += rate * getFactor<s, Tup>())... };
	(void)x;
	//return res;
}

template <typename Tup>
void increaseProducts(Tup &tup, float rate)
{
	increaseProducts(tup, rate, std::make_index_sequence<std::tuple_size<Tup>::value>());
}

template <typename Tup, size_t a, size_t...s>
typename std::enable_if<(std::index_sequence<s...>::size() == 0u), float>::type minConcentration(Tup const &tup, std::index_sequence<a>)
{
	return std::get<a>(tup).destination;
}

template <typename Tup, size_t a, size_t...s>
typename std::enable_if<(std::index_sequence<s...>::size() > 0u), float>::type minConcentration(Tup const &tup, std::index_sequence<a, s...>)
{
	auto res = std::get<a>(tup).destination;
	auto x = { (res < std::get<s>(tup).destination ? true : res = std::get<s>(tup).destination)... };
	(void)x;
	return res;
}

template <typename Tup>
float findMinConcentration(Tup const &tup)
{
	return minConcentration(tup, std::make_index_sequence<std::tuple_size<Tup>::value>());
}

template <typename Tup, size_t a, size_t...s>
constexpr typename std::enable_if<(std::index_sequence<s...>::size() == 0u), float>::type maxFactor(const Tup tup, const std::index_sequence<a>)
{
	return getFactor<a, Tup>();
}

template <typename Tup, size_t a, size_t...s>
constexpr typename std::enable_if<(std::index_sequence<s...>::size() > 0u), float>::type maxFactor(const Tup tup, const std::index_sequence<a, s...>)
{
	return getFactor<a, Tup>() > maxFactor(tup, std::index_sequence<s...>()) ? getFactor<a, Tup>() : maxFactor(tup, std::index_sequence<s...>());
}

template <typename Tup>
constexpr float findMaxFactor(const Tup tup)
{
	return maxFactor(tup, std::make_index_sequence<std::tuple_size<Tup>::value>());
}

template<typename substrate_tuple, typename product_tuple>
inline void productInhibitionEquation(float rate, substrate_tuple st, product_tuple pt)
{
	// First calculate actual rate
	auto substrate_val = findMinConcentration(st);	
	auto dividend = rate * substrate_val;
  auto sum_val = sumTuple(pt);
	auto divisor = rate + substrate_val + sum_val;
  auto factor = findMaxFactor(st);
	auto change = dividend / (divisor * factor);
	// Second modify values
	decreaseSubstrates(st, change );
	increaseProducts(pt, change );
}
template<typename substrate_tuple, typename product_tuple>
inline void simpleEquation(float rate, substrate_tuple st, product_tuple pt)
{
	// First calculate actual rate
	auto substrate_val = findMinConcentration(st);
	auto dividend = rate * substrate_val;
	auto divisor = rate + substrate_val;
  auto factor = findMaxFactor(st);
	auto change = dividend / (divisor * factor);
	// Second modify values
	decreaseSubstrates(st, change);
	increaseProducts(pt, change);
}
template<typename sub_1, typename sub_2>
inline void clampDefusion(const float rate, sub_1 sub1, sub_2 sub2, float volume_1, float volume_2)
{
	auto sub_1_conc = sub1.destination / volume_1;
	auto sub_2_conc = sub2.destination / volume_2;
	auto ratio = sub_1::factor / sub_2::factor;
	auto diff = sub_2_conc - sub_1_conc;
	auto change = rate * diff * volume_1;
	sub1.destination += change;
	sub2.destination -= change * ratio;
}
template<typename sub_1, typename sub_2>
inline void logDefusion(const float rate, sub_1 sub1, sub_2 sub2, float volume_1, float volume_2)
{
	auto sub_1_conc = sub1.destination / volume_1;
	auto sub_2_conc = sub2.destination / volume_2;
	auto diff = sub_2_conc - sub_1_conc;
	auto dividend = rate * diff;
	auto divisor = rate + std::abs(diff);
	auto change = dividend / divisor * volume_1;
	auto ratio = sub_1::factor / sub_2::factor;
	sub1.destination += change;
	sub2.destination -= change * ratio;
}
}

