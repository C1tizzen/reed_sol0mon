#include <cmath>
#include <vector>
#include <iostream>
#include <stdint.h>
#include <string.h>
#include "gf.h"


uint8_t exp_h[512]; // Для эффективности gf_exp [] имеет размер 2 * GF_SIZE, так что простое умножение двух чисел может быть решено без вызова% 255
uint8_t log_h[256];

#define GENERATOR 2 // заменяем альфу в классической теории

// Здесь мы используем алгоритм умножения по методу «русских мужиков» для нахождения неуменьшенного многочлена
/*@ краткое умножение в полях Галуа
* @param x - левый операнд
* @param y - правильный операнд
* @param prim - примитивный двоичный многочлен
* @param field_charac_full - количество элементов в поле Галуа
* @return x * y */
int gf_mult(int x, int y, int prim = 0, int field_charac_full=256) {
	int res = 0;

	while (y) {
		if ((y & 1) == 1) {
			res = res ^ x;
		}
		y = y >> 1; //~ y/2
		x = x << 1; //~ x*2

		if ((prim > 0) && (x & field_charac_full)) { // если x ~ 256 использовать xor
			x = x ^ prim;
		}
		
	}

	return res;
}


/*@brief Нахождение списка всех простых многочленов для заданного генератора(2) и характеристического показателя поля Галуа(8)
* @param c_exp - показатель степени поля Галуа
* @return[неприводимый многочлен] */
std::vector<int> find_prime_polys(int c_exp = 8) {
	int root_charac = 2;
	int field_charac = pow(root_charac, c_exp) - 1;             //255
	int field_charac_next = pow(root_charac, c_exp + 1) - 1;    // 511

	std::vector<int> prim_candidates;
	std::vector<int> irreducible_polynomial;
	std::vector<int> seen;

	for (int j = 0; j < field_charac + 1; j++) {
		seen.push_back(0);
	}

	for (int i = field_charac + 2; i < field_charac_next; i += root_charac) { // нам нужно число> 256, чтобы не повторяться
		prim_candidates.push_back(i);                                         // пробуем каждый простой многочлен, исключая четные
	}

// Здесь реализован метод грубой силы, чтобы найти все эти простые многочлены, генерируя все возможные
// простые многочлены (т.е. все целые числа между field_charac + 1 и field_charac * 2), а затем мы строим
// все Поле Галуа, и мы отклоняем кандидатный простой многочлен, если он дублирует хотя бы одно значение или если он
// генерирует значение выше field_charac (т. е. вызывает переполнение).
	for (int i = 0; i < prim_candidates.size(); i++) {
		int prim = prim_candidates[i];
		for (int j = 0; j < field_charac + 1; j++) {
			seen[j] = 0;                                       // указывает, было ли значение в поле уже сгенерировано (замечено [x]? == 1) или нет
		}

		bool conflict = false;                                 // чтобы узнать, был ли хотя бы 1 конфликт
		int x = 1;

		for (int j = 0; j < field_charac; j++) {
			x = gf_mult(x, GENERATOR, prim, field_charac + 1); // вычисляем следующее значение в поле

			if (x > field_charac || seen[x] == 1) {            // если это число - дубликат, то мы его отклоняем
				//std::cout << x << std::endl;
				conflict = true;
				break;
			}
			else {
				seen[x] = 1;                                   // запоминаем это значение для обнаружения будущих дубликатов
			}
		}

		if (!conflict) {                                       // если конфликта не было, то это простой многочлен
			irreducible_polynomial.push_back(prim);
		}

	}

	return irreducible_polynomial;
}

/*@brief Вычисление таблиц "логарифм" и "экспоненты" для дальнейшего использования в арифметических операциях
* @param gf_exp - переданный по ссылке вектор, в котором будут помещены значения exp
* @param gf_log - переданный по ссылке вектор, в котором будут помещены значения журнала
* @param prim - примитивный двоичный многочлен
* @param c_exp - показатель степени поля Галуа */
void init_tables(std::vector<int>& gf_exp, std::vector<int>& gf_log, int prim = 285, int c_exp = 8) {
	int field_charac = pow(2, c_exp) - 1;        //255

	for (int i = 0; i < field_charac * 2; i++) { //инизиализация
		gf_exp.push_back(0);
	}
	for (int i = 0; i < field_charac + 1; i++) {
		gf_log.push_back(0);
	}

	// Для каждого элемента из поля Галуа мы вычисляем log и exp
	int x = 1;
	for (int i = 0; i < field_charac; i++) {
		gf_exp[i] = x;
		gf_log[x] = i;
		x = gf_mult(x, GENERATOR, prim, field_charac + 1);
	}

	for (int i = field_charac; i < field_charac * 2; i++) {
		gf_exp[i] = gf_exp[i - field_charac];
	}
}



