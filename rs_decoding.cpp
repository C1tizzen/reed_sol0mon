#include <iostream>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include "gf.h"

#define GENERATOR 2 // заменяем альфу в классической теории


/*@brief Функция вычисляет полином синдромов ошибок
* @param msg - входящее сообщение, которое представлено вектором многочленов(в данном случае целых чисел)
* @param red_code_len - количество символов, представляющих избыточный код
* @ полином генератора возврата*/
std::vector<int> rs_calc_syndromes(std::vector<int>& msg, int& red_code_len) {
	std::vector<int> synd;
	synd.reserve(red_code_len + 1);
	
	for (int i = 0; i < red_code_len + 1; i++) {     //Инициализация
		synd.push_back(0);
	}

	for (int i = 1; i < red_code_len + 1; i++) {    // Нам нужен первый элемент 0 для математической точности, иначе будут ошибки
		int temp = gf::pow_my(GENERATOR, i - 1);
		synd[i] = gf::gf_poly_eval(msg, temp);       //S[i+1] = C(a^(i)) 
	}

	return synd;
}

/*@brief Функция найти полином локатора ошибок L(x) в соответствии с известным местом ошибки(для исправления ошибки)
* @param err_pos - вектор с позициями ошибочных символов
* @ возвратный многочлен локатора ошибок L(x) */
std::vector<int> rs_find_errarta_locator(std::vector<int>& err_pos) {
	std::vector<int> e_loc, temp, temp2, add;
	e_loc.reserve(err_pos.size());
	temp.reserve(err_pos.size());
	temp2.reserve(err_pos.size());
	add.reserve(err_pos.size());

	e_loc.push_back(1);

	temp.push_back(1);

	temp2.push_back(0);
	temp2.push_back(0);

	for (int i = 0; i < err_pos.size(); i++) {
		temp2[0] = gf::pow_my(GENERATOR, err_pos[i]);         // поскольку мы знаем местоположение ошибки, мы можем найти L (x) как
		add = gf::gf_poly_add(temp, temp2);                   // L(x) = П(1 + x*alpha^(i))
		e_loc = gf::gf_poly_mult(e_loc, add);
	}

	return e_loc;
}


/*@brief Функция найти полином ошибок W(x) = S(x) * L(x)
* @param Syndrome - многочлен синдромов ошибок(вектор int) S(x)
* @param err_loc - многочлен локатора ошибок L(x)
* @param err_loc_size - размер L(x)
* @ полином ошибки возврата W(x) */
std::vector<int> rs_find_error_evaluator(std::vector<int> synd, std::vector<int> err_loc, int err_loc_size) {
	std::vector<int> poly_mul, remainder;
	poly_mul.reserve(synd.size() + err_loc.size());

	poly_mul = gf::gf_poly_mult(synd, err_loc);
	remainder.reserve(err_loc_size);
	
	for (int i = poly_mul.size() - err_loc_size; i < poly_mul.size(); i++) {
		remainder.push_back(poly_mul[i]);
	}

// Поскольку W (x) не может превышать u-1, где u - количество ошибок, мы будем использовать деление, чтобы отбросить лишнюю часть
// gf :: gf_poly_div (poly_mul, справка, частное, остаток);

	return remainder;
}


/*@brief Исправление входящего сообщения с использованием алгоритма Форни, который вычисляет значение ошибки
* @param msg_in - входящее сообщение, которое представлено вектором многочленов(в данном случае целых чисел)
* @param Syndrome - многочлен синдромов ошибок(вектор int) S(x)
* @param err_pos - вектор с позициями ошибочных символов
* @return исправленное сообщение ввода */
std::vector<int> rs_correct_errata(std::vector<int> msg_in, std::vector<int> synd, std::vector<int> err_pos) {
	std::vector<int> coef_pos, err_loc, err_eval, x;
	coef_pos.reserve(err_pos.size());
	err_loc.reserve(err_pos.size());
	err_eval.reserve(err_pos.size());
	x.reserve(coef_pos.size());

	int len = msg_in.size();
	for (int i = 0; i < err_pos.size(); i++) {                         // конвертируем ошибочные позиции в коэффициенты степени
		coef_pos.push_back(len - 1 - err_pos[i]);
	}


	// находим полином локатора ошибок L (x) в соответствии с известным местоположением ошибки
	err_loc = rs_find_errarta_locator(coef_pos);
	
	// находим полином ошибки WITH (x)
	reverse(synd.begin(), synd.end());
	err_eval = rs_find_error_evaluator(synd, err_loc, err_loc.size());
	reverse(err_eval.begin(), err_eval.end());

	// x - сохранит позицию ошибок
    // нам нужно получить полином X местоположения ошибки из позиций ошибки в err_pos
    // (корни многочлена локатора ошибок, т. е. где он принимает значение 0)
	for (int i = 0; i < coef_pos.size(); i++) {
		int l = 255 - coef_pos[i];
		x.push_back(gf::pow_my(GENERATOR, -l));
	}

	// с помощью алгоритма Форни находим значения ошибок
	//Алгоритм Форни используется для оценки значений ошибок. 
	//Здесь положение ошибки и коэффициенты ω(x) берутся здесь как вход. 
	//Он также использует множитель поля Галуа в качестве алгоритма поиска Ченя.
	
	
	std::vector<int> E, err_loc_prime_tmp;
	E.reserve(msg_in.size());
	err_loc_prime_tmp.reserve(err_pos.size());

	for (int i = 0; i < msg_in.size(); i++) {
		E.push_back(0);                                 // сохраним значения, которые необходимо исправить, в исходное сообщение с ошибками
	}


	int Xlength = x.size();
	reverse(err_eval.begin(), err_eval.end());

	for (int i = 0; i < x.size(); i++) {
		int Xi_inv = gf::inverse(x[i]);

		// Находим формальную производную многочлена локатора ошибок
// формальная производная локатора ошибок используется в качестве знаменателя алгоритма Форни,
// который просто говорит, что i-е значение ошибки задается error_evaluator (gf_inverse (Xi)) / error_locator_derivative (gf_inverse (Xi)).
		for (int j = 0; j < Xlength; j++) {
			if (j != i) {
				err_loc_prime_tmp.push_back(gf::sub(1, gf::mult(Xi_inv, x[j])));
			}
		}

		// многочлен ошибок Yi = W (Xi ^ (- 1)) / L '(Xi ^ (- 1))
// вычисляем произведение, которое является знаменателем алгоритма Форни (производная локатора ошибок)
		int err_loc_prime = 1;

		for (int coef = 0; coef < err_loc_prime_tmp.size(); coef++) {
			err_loc_prime = gf::mult(err_loc_prime, err_loc_prime_tmp[coef]);
		}
		
		err_loc_prime_tmp.clear();

		int y;                                              
		y = gf::gf_poly_eval(err_eval, Xi_inv);                          // числитель
		y = gf::mult(gf::pow_my(x[i], 1), y);

		try {
			if (err_loc_prime == 0) {                                     // делитель не должен быть 0
				throw "Could not find error magnitude";
			}
		}
		catch (const char* ex) {
			std::cerr << std::endl << std::endl << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =";
			std::cerr << std::endl << "Error: " << ex << std::endl;
			std::cerr << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" << std::endl;
			exit(-1);
		}

		int magnitude = gf::div_my(y, err_loc_prime);                      //The value of the error value calculated by the Forney algorithm
		                                                                   //Dividing the error estimator by the derivative of the error locator 
		E[err_pos[i]] = magnitude;                                         //gives us the error value, that is, the value for recovering the symbol
	}

	// исправляем наш ошибочный многочлен
	msg_in = gf::gf_poly_add(msg_in, E);                                   // C(x) = C'(x) + E(x) (xor)
	return msg_in;
}


/*@brief Функция найти полином локатора ошибок L(x) с использованием алгоритма Берлекампа - Месси
* @param Syndrome - многочлен синдромов ошибок(вектор int) S(x)
* @param red_code_len - количество символов, представляющих избыточный код
* @ возвратный многочлен локатора ошибок L(x)*/
std::vector<int> rs_find_error_locator(std::vector<int> synd, int red_code_len) {
	std::vector<int> err_loc, old_loc;
	err_loc.reserve(synd.size());
	old_loc.reserve(synd.size());

	err_loc.push_back(1);                               // многочлен локатора ошибок(sigma) C(x)
	old_loc.push_back(1);                               // полином локатора ошибок предыдущей итерации
	
	
	int synd_shift = synd.size() - red_code_len;

// Алгоритм Берлекампа – Месси является альтернативой декодеру Рида – Соломона Петерсона для решения системы линейных уравнений.
// Основная идея заключается в том, что алгоритм итеративно вычисляет полином локатора ошибок. Для этого он вычисляет дельта-расхождение,
// по которому мы можем определить, нужно ли нам обновлять локатор ошибок или нет
	int k = 0, delta = 0;
	for (int i = 0; i < red_code_len; i++) {
		k = i + synd_shift;
		
		// вычисляем дельту несоответствия
		delta = synd[k];
				
		for (int j = 1; j < err_loc.size() ; j++) {
			delta ^= gf::mult(err_loc[err_loc.size() - 1 - j], synd[k - j]);    // delta = Sn + C1*Sn-1 +..+ Cj*Sk-j
		}

		// сдвигаем многочлены для вычисления следующей степени
		old_loc.push_back(0);
		
		std::vector<int> new_loc;
		if(delta != 0){                                                         // если delta == 0, алгоритм считает, что C (x) и L верны на данный момент, и продолжает
			if (old_loc.size() > err_loc.size()) {                              //~2*L <= k + erase_count
				// Вычисление полинома локатора ошибок Sigma
				new_loc = gf::gf_poly_scale(old_loc, delta);
				old_loc = gf::gf_poly_scale(err_loc, gf::inverse(delta));
				err_loc = new_loc;
			}

			// Обновление с расхождением
			err_loc = gf::gf_poly_add(err_loc, gf::gf_poly_scale(old_loc, delta));

		}
	}

	while (err_loc.size() && err_loc[0] == 0) {
		err_loc.erase(err_loc.begin());
	}


	int errs = err_loc.size() - 1;
	try {
		if (errs*2 > red_code_len) {
			throw "Too many errors to correct";
		}
	}
	catch (const char* ex) {
		std::cerr << std::endl << std::endl << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =";
		std::cerr << std::endl << "Error: " << ex << std::endl;
		std::cerr << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" << std::endl;
		exit(-1);
	}

	return err_loc;
}


/*@brief Функция находим нули этого многочлена
* @param err_loc - многочлен локатора ошибок L(x)
* @param nmess - размер сообщения
* @return вектор индекса символов, которые необходимо исправить*/
std::vector<int> rs_find_errors(std::vector<int> err_loc, int nmess) {
// стираем с помощью полинома локатора ошибок, мы просто используем пробную замену, чтобы найти нули этого полинома,
// который определяет места ошибки (т. е. индекс символов, которые необходимо исправить)
	std::vector<int> err_pos;
	err_pos.reserve(nmess);

	int errs = err_loc.size() - 1;
	for (int i = 0; i < nmess; i++) {
		if (gf::gf_poly_eval(err_loc, gf::pow_my(GENERATOR, i)) == 0) {  // если 0, то это корень многочлена локатора ошибок
			err_pos.push_back(nmess - 1 - i);
		}
	}

	try {
		if (err_pos.size() != errs) {
			throw "Too many (or few) errors found for the errata locator polynomial!";
		}
	}
	catch (const char* ex) {
		std::cerr << std::endl << std::endl << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =";
		std::cerr << std::endl << "Error: " << ex << std::endl;
		std::cerr << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" << std::endl;
		exit(-1);
	}
	
	return err_pos;
}


/*@brief Функция декодирования сообщения
* @param msg_in - вектор закодированного сообщения
* @param red_code_len - количество символов, представляющих избыточный код
* @param erase_pos - известные ошибки позиции
* @return декодированное сообщение */
std::vector<int> rs_decode_msg(std::vector<int> msg_in, int red_code_len) {
	try {
		if (msg_in.size() > 255) {
			throw "Message is too long";
		}
	}
	catch (const char* ex) {
		std::cerr << std::endl << std::endl << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =";
		std::cerr << std::endl << "Error: " << ex << std::endl;
		std::cerr << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" << std::endl;
		exit(-1);
	}

	std::vector<int> msg_out = msg_in;

	// чтобы мы не считали порождающий многочлен несколько раз и не делили,
	// сразу подсчитываем полином синдрома ошибки, и если его нет хотя бы
	// в нем одно ненулевое значение, тогда сообщение не искажается
	std::vector<int> synd = rs_calc_syndromes(msg_out, red_code_len);
	int max = *max_element(synd.begin(), synd.end());

	if (max == 0) {
		return msg_out;
	}


	std::vector<int> err_pos, err_loc;
	err_loc.reserve(synd.size());
	err_pos.reserve(synd.size());

	// Находим полином локатора ошибок L (x)
		err_loc = rs_find_error_locator(synd, red_code_len);
		
		reverse(err_loc.begin(), err_loc.end());

		// находим вектор индекса символов, которые нужно исправить
		err_pos = rs_find_errors(err_loc, msg_out.size());
		
		try {
			if (err_pos.size() == 0) {
				throw "Could not locate error";
			}
		}
		catch (const char* ex) {
			std::cerr << std::endl << std::endl << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =";
			std::cerr << std::endl << "Error: " << ex << std::endl;
			std::cerr << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" << std::endl;
			exit(-1);
		}

	
		// Находим значения ошибок и применяем их для исправления сообщения
	// вычисление полиномов оценки и количества ошибок, затем исправление ошибок и стирания
	msg_out = rs_correct_errata(msg_out, synd, err_pos);
	
	// подсчитываем полином синдрома ошибки, и если его нет хотя бы
// в нем одно ненулевое значение, тогда сообщение декодируется успешно
	synd = rs_calc_syndromes(msg_out, red_code_len);
	max = *max_element(synd.begin(), synd.end());
	try {
		if (max > 0) {
			throw "Could not correct message";
		}
	}
	catch (const char* ex) {
		std::cerr << std::endl << std::endl << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =";
		std::cerr << std::endl << "Error: " << ex << std::endl;
		std::cerr << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" << std::endl;
		exit(-1);
	}

	return msg_out;
}
