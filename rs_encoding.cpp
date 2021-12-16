#include <iostream>
#include <string.h>
#include <stdint.h>
#include <vector>
#include "gf.h"


#define GENERATOR 2 // �������� ����� � ������������ ������

// �������, ������� ��������� ������������ ������� ��� ��������� ���������� ������������� ��������
/* @brief ��������� ����������� �������
* @param red_code_len - ���������� ��������, �������������� ���������� ���
* @ ������� ���������� �������� */
	std::vector<int> rs_generator_poly(int &red_code_len) {
		std::vector<int> g, temp;
		g.reserve(red_code_len);
		temp.reserve(red_code_len);

		g.push_back(1);
		temp.push_back(1);
		temp.push_back(0);

		for (uint8_t i = 0; i < red_code_len; i++) {
			temp[1] = gf::pow_my(GENERATOR, i);
			g = gf::gf_poly_mult(g, temp);
		}
		return g;                                         //(x-1)*(x-2^(1))*..*(x-2^(red_code_len-1))
	}


	/* @ ������� ��������� ���������
		* @param msg_in - �������� ���������, �������������� �������� �����������(� ������ ������ ����� �����)
		* @param red_code_len - ���������� ��������, �������������� ���������� ���
		* @return �������������� ��������� = ������[msg_in] + [���������� ����������](� ������ ������ ����� �����) */
	std::vector<int> rs_encode_msg(std::vector<int>& msg_in, int red_code_len) {
		try {
			if (msg_in.size() + red_code_len < 256) {
				int msg_in_size = msg_in.size();
				std::vector<int> gen;
				gen.reserve(red_code_len);

				gen = rs_generator_poly(red_code_len);

				std::vector<int> msg_out;
				int msg_out_size = msg_in_size + red_code_len;
				msg_out.reserve(msg_out_size);

				for (int i = 0; i < msg_out_size; i++) {// �������������� msg_out: len = msg_in.size () + red_code_len
					if (i < msg_in_size)
						msg_out.push_back(msg_in[i]);                      // ������ k �������� �������� �������� ���������
					else
					{
						msg_out.push_back(0);
					}
				}

				std::vector<int> quotient, remainder;
				quotient.reserve(msg_out_size);
				remainder.reserve(msg_out_size);

				gf::gf_poly_div(msg_out, gen, quotient, remainder);        // ����� �������� ��������� �� ����������� � ���������� �������
				reverse(remainder.begin(), remainder.end());

				for (int i = 0; i < msg_out_size; i++) { 
					if (i < msg_in_size)
						msg_out[i] = msg_in[i];                      
					else
					{
						msg_out[i] = remainder[msg_out_size - 1 - i];
					}
				}

// � ������� ������� ������� ��� ����������� �����, ��� ��� ��� ���� ������ ������ ���������, �� ��� ���� ��������!

//// ����� �������������� �������
// for (int i = 0; i <msg_in.size (); i ++) {
// int coef = msg_out [i];
// if (coef! = 0) {// log [0] �� ����������, ������� ��� ����� ������� ��������� ���� ������
// for (int j = 1; j <gen.size (); j ++) {
// msg_out [i + j] ^ = gf :: mul (gen [j], coef); // ~ '+ ='
//}
//}
//}

// for (int i = 0; i <msg_in.size (); i ++) {// � ������� k �������� ������ ���� �������,
// msg_out [i] = msg_in [i]; // �� ��� �� ����� ��� �����������, ������� ��������� ����� msg_in
//}

				return msg_out;
			}
			else {
				throw "\"The total number of characters - messages + redundant code - exceeds 256\"";
			}
		}
		catch (const char* ex) {
			std::cerr << std::endl << std::endl << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =";
			std::cerr << std::endl << "Error: " << ex << std::endl;
			std::cerr << "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" << std::endl;
		}
	}