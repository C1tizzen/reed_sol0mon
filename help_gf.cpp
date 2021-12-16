#include <cmath>
#include <vector>
#include <iostream>
#include <stdint.h>
#include <string.h>
#include "gf.h"


uint8_t exp_h[512]; // ��� ������������� gf_exp [] ����� ������ 2 * GF_SIZE, ��� ��� ������� ��������� ���� ����� ����� ���� ������ ��� ������% 255
uint8_t log_h[256];

#define GENERATOR 2 // �������� ����� � ������������ ������

// ����� �� ���������� �������� ��������� �� ������ �������� ������� ��� ���������� �������������� ����������
/*@ ������� ��������� � ����� �����
* @param x - ����� �������
* @param y - ���������� �������
* @param prim - ����������� �������� ���������
* @param field_charac_full - ���������� ��������� � ���� �����
* @return x * y */
int gf_mult(int x, int y, int prim = 0, int field_charac_full=256) {
	int res = 0;

	while (y) {
		if ((y & 1) == 1) {
			res = res ^ x;
		}
		y = y >> 1; //~ y/2
		x = x << 1; //~ x*2

		if ((prim > 0) && (x & field_charac_full)) { // ���� x ~ 256 ������������ xor
			x = x ^ prim;
		}
		
	}

	return res;
}


/*@brief ���������� ������ ���� ������� ����������� ��� ��������� ����������(2) � ������������������� ���������� ���� �����(8)
* @param c_exp - ���������� ������� ���� �����
* @return[������������ ���������] */
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

	for (int i = field_charac + 2; i < field_charac_next; i += root_charac) { // ��� ����� �����> 256, ����� �� �����������
		prim_candidates.push_back(i);                                         // ������� ������ ������� ���������, �������� ������
	}

// ����� ���������� ����� ������ ����, ����� ����� ��� ��� ������� ����������, ��������� ��� ���������
// ������� ���������� (�.�. ��� ����� ����� ����� field_charac + 1 � field_charac * 2), � ����� �� ������
// ��� ���� �����, � �� ��������� ����������� ������� ���������, ���� �� ��������� ���� �� ���� �������� ��� ���� ��
// ���������� �������� ���� field_charac (�. �. �������� ������������).
	for (int i = 0; i < prim_candidates.size(); i++) {
		int prim = prim_candidates[i];
		for (int j = 0; j < field_charac + 1; j++) {
			seen[j] = 0;                                       // ���������, ���� �� �������� � ���� ��� ������������� (�������� [x]? == 1) ��� ���
		}

		bool conflict = false;                                 // ����� ������, ��� �� ���� �� 1 ��������
		int x = 1;

		for (int j = 0; j < field_charac; j++) {
			x = gf_mult(x, GENERATOR, prim, field_charac + 1); // ��������� ��������� �������� � ����

			if (x > field_charac || seen[x] == 1) {            // ���� ��� ����� - ��������, �� �� ��� ���������
				//std::cout << x << std::endl;
				conflict = true;
				break;
			}
			else {
				seen[x] = 1;                                   // ���������� ��� �������� ��� ����������� ������� ����������
			}
		}

		if (!conflict) {                                       // ���� ��������� �� ����, �� ��� ������� ���������
			irreducible_polynomial.push_back(prim);
		}

	}

	return irreducible_polynomial;
}

/*@brief ���������� ������ "��������" � "����������" ��� ����������� ������������� � �������������� ���������
* @param gf_exp - ���������� �� ������ ������, � ������� ����� �������� �������� exp
* @param gf_log - ���������� �� ������ ������, � ������� ����� �������� �������� �������
* @param prim - ����������� �������� ���������
* @param c_exp - ���������� ������� ���� ����� */
void init_tables(std::vector<int>& gf_exp, std::vector<int>& gf_log, int prim = 285, int c_exp = 8) {
	int field_charac = pow(2, c_exp) - 1;        //255

	for (int i = 0; i < field_charac * 2; i++) { //�������������
		gf_exp.push_back(0);
	}
	for (int i = 0; i < field_charac + 1; i++) {
		gf_log.push_back(0);
	}

	// ��� ������� �������� �� ���� ����� �� ��������� log � exp
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



