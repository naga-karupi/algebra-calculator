/**
 * @file multi-mat.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2025-04-17
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#include <iostream>
#include <cstdlib>

class Matrix {
	private:
		int rows;
		int columns;
		int *v;

	public:
		Matrix(int sizeRows, int sizeColmns, int *data);
		Matrix(int sizeRows, int sizeColmns);
		void freeMat();
		int getMat(int x,  int y);
		void setMat(int x, int y, int val);

		void printMat();

		void Strassen(Matrix *matA, Matrix *matB, Matrix *matC);

		static void errorMat(const char *str);
};

Matrix::Matrix(int sizeRows, int sizeColmns, int *data) {
	rows = sizeRows;
	columns = sizeColmns;
	v = data;
}

Matrix::Matrix(int sizeRows, int sizeColmns) {
	rows = sizeRows;
	columns = sizeColmns;
}

void Matrix::freeMat() {
	free(v);
	free(this);
}

int Matrix::getMat(int x, int y) {
	if(x >= rows) errorMat("getMat: x is out of range");
	if(y >= columns) errorMat("getMat: y is out of range");
	return v[x*this->columns + y];
}

void Matrix::setMat(int x, int y, int val) {
	if(x >= rows) errorMat("setMat: x is out of range");
	if(y >= columns) errorMat("setMat: y is out of range");
	v[x*this->columns + y] = val;
}

void Matrix::printMat(){
	for(int i=0; i<rows; i++) {
		for(int j=0; j<columns; j++) {
			std::cout << " " << getMat(i, j);
		}
		std::cout << std::endl;
	}
}

void Matrix::Strassen(Matrix *matA, Matrix *matB, Matrix *matC) {
	if(matA->rows == 1) {
		for(int i=0; i<matB->columns; i++) {
			for(int j=0;j<matA->columns; j++) {
				matC->v[i] += matC->v[i] + matA->v[j] * matB->v[j+i*matB->rows];
			}
		}
		return;
	}
	else if(matB->columns == 1) {
		for(int i=0; i<matA->rows; i++) {
			for(int j=0; j<matB->rows; j++) {
				matC->v[i] += matA->v[j+i*matA->columns] * matB->v[j];
			}
		}
		return;
	}
	else if(matA->columns == 1 && matB->rows == 1) {
		for(int i=0; i<matA->rows; i++) {
			for(int j=0; j<matB->columns; j++) {
				matC->v[j+i*matB->columns] = matA->v[j]*matB->v[i];
			}
		}
		return;
	}
	
	Matrix *matA11, *matA12, *matA21, *matA22;
	Matrix *matB11, *matB12, *matB21, *matB22;
	Matrix *matC11, *matC12, *matC21, *matC22;
	
	int i = 0, j = 0, x, y;
	int sizeARow1, sizeARow2, sizeAColumn1, sizeAColumn2;
	int sizeBRow1, sizeBRow2, sizeBColumn1, sizeBColumn2;
	sizeARow1 	 = matA->rows / 2;
	sizeARow2 	 = matA->rows - sizeARow1;
	sizeAColumn1 = matA->columns / 2;
	sizeAColumn2 = matA->columns - sizeAColumn1;
	sizeBRow1	 = sizeAColumn1;
	sizeBRow2 	 = sizeAColumn2;
	sizeBColumn1 = matB->columns / 2;
	sizeBColumn2 = matB->columns - sizeBColumn1;
	
	matA11 = new Matrix(sizeARow2, sizeAColumn2); matA12 = new Matrix(sizeARow2, sizeAColumn2);
	matA21 = new Matrix(sizeARow2, sizeAColumn2); matA22 = new Matrix(sizeARow2, sizeAColumn2);
	matB11 = new Matrix(sizeBRow2, sizeBColumn2); matB12 = new Matrix(sizeBRow2, sizeBColumn2);
	matB21 = new Matrix(sizeBRow2, sizeBColumn2); matB22 = new Matrix(sizeBRow2, sizeBColumn2);
	matC11 = new Matrix(sizeARow2, sizeBColumn2); matC12 = new Matrix(sizeARow2, sizeBColumn2);
	matC21 = new Matrix(sizeARow2, sizeBColumn2); matC22 = new Matrix(sizeARow2, sizeBColumn2);


	int val;
	for(x=0; x<sizeARow1; x++) {
		for(y=0; y<sizeAColumn1; y++) {
			val = matA->getMat(x, y);
			matA11->setMat(x, y, val);
		}
	}
	
	for(x=0; x<sizeARow1; x++) {
		for(y=sizeAColumn1; y<matA->columns; y++) {
			val = matA->getMat(x, y);
			matA12->setMat(x, j, val);
			j++;
		}
		j = 0;
	}

	for(x=sizeARow1; x<matA->rows; x++) {
		for(y=0; y<sizeAColumn1; y++) {
			val = matA->getMat(x, y);
			matA21->setMat(i, y, val);
		}
		i++;
	}
	i = 0; j = 0;

	for(x=sizeARow1; x<matA->rows; x++) {
		for(y=sizeAColumn1; y<matA->columns; y++) {
			val = matA->getMat(x, y);
			matA22->setMat(i, j, val);
			j++;
		}
		j = 0;
		i++;
	}
	i = 0;

	for(x=0; x<sizeBRow1; x++) {
		for(y=0; y<sizeBColumn1; y++) {
			val = matB->getMat(x, y);
			matB11->setMat(x, y, val);
		}
	}

	for(x=0; x<sizeBRow1; x++) {
		for(y=sizeBColumn1; y<matB->columns; y++) {
			val = matB->getMat(x, y);
			matB12->setMat(x, j, val);
			j++;
		}
		j = 0;
	}

	for(x=sizeBRow1; x<matB->rows; x++) {
		for(y=0; y<sizeBColumn1; y++) {
			val = matB->getMat(x, y);
			matB21->setMat(i, y, val);
		}
		i++;
	}
	i = 0;
	
	for(x=sizeBRow1; x<matB->rows; x++) {
		for(y=sizeBColumn1; y<matB->columns; y++) {
			val = matB->getMat(x, y);
			matB22->setMat(i, j, val);
			j++;
		}
		i++;
		j = 0;
	}
	
	//conquer
	int val1, val2, val3, val4, val5, val6, val7;
   	Matrix *matP1, *matP2, *matP3, *matP4, *matP5, *matP6, *matP7;
   	Matrix *sum1, *sum2, *sub1, *sub2;
	// エラー
	matP1 = new Matrix(sizeARow2, sizeBColumn2);
	matP2 = new Matrix(sizeARow2, sizeBColumn2);
	matP3 = new Matrix(sizeARow2, sizeBColumn2);
	matP4 = new Matrix(sizeARow2, sizeBColumn2);
	matP5 = new Matrix(sizeARow2, sizeBColumn2);
	matP6 = new Matrix(sizeARow2, sizeBColumn2);
	matP7 = new Matrix(sizeARow2, sizeBColumn2);
	sum1 = new Matrix(sizeARow2, sizeBColumn2);
	sum2 = new Matrix(sizeARow2, sizeBColumn2);
	sub1 = new Matrix(sizeARow2, sizeBColumn2);
	sub2 = new Matrix(sizeARow2, sizeBColumn2);
	for(x=0; x<sizeARow2; x++) {
		for(y=0; y<sizeAColumn2; y++) {
			val1 = matA11->getMat(x, y);
			val2 = matA22->getMat(x, y);
			sum1->setMat(x, y, val1+val2);
		}
	}

	for(x=0; x<sizeBRow1; x++) {
		for(y=0; y<sizeBColumn2; y++) {
			val3 = matB11->getMat(x, y);
			val4 = matB22->getMat(x, y);
			sum2->setMat(x, y, val3+val4);
		}
	}
	Strassen(sum1, sum2, matP1);

	for(x=0; x<sizeARow2; x++) {
		for(y=0; y<sizeAColumn2; y++) {
			val1 = matA21->getMat(x, y);
			val2 = matA22->getMat(x, y);
			sum1->setMat(x, y, val1+val2);
		}
	}
	Strassen(sum1, matB11, matP2);

	for(x=0; x<sizeBRow2; x++) {
		for(y=0; y<sizeBColumn2; y++) {
			val3 = matB12->getMat(x, y);
			val4 = matB22->getMat(x, y);
			sub2->setMat(x, y, val3-val4);
		}
	}
	Strassen(matA11, sub2, matP3);

	for(x=0; x<sizeBRow2; x++) {
		for(y=0; y<sizeBColumn2; y++) {
			val3 = matB21->getMat(x, y);
			val4 = matB11->getMat(x, y);
			sub2->setMat(x, y, val3-val4);
		}
	}
	Strassen(matA22, sub2, matP4);

	for(x=0; x<sizeBRow2; x++) {
		for(y=0; y<sizeBColumn2; y++) {
			val1 = matA11->getMat(x, y);
			val2 = matA12->getMat(x, y);
			sum1->setMat(x, y, val1+val2);
		}
	}
	Strassen(sum1, matB22, matP5);

	for(x=0; x<sizeBRow2; x++) {
		for(y=0; y<sizeBColumn2; y++) {
			val1 = matA21->getMat(x, y);
			val2 = matA11->getMat(x, y);
			sub1->setMat(x, y, val1-val2);
		}
	}

	for(x=0; x<sizeBRow2; x++) {
		for(y=0; y<sizeBColumn2; y++) {
			val3 = matB11->getMat(x, y);
			val4 = matB12->getMat(x, y);
			sum2->setMat(x, y, val3+val4);
		}
	}
	Strassen(sub1, sum2, matP6);

	for(x=0; x<sizeBRow2; x++) {
		for(y=0; y<sizeBColumn2; y++) {
			val1 = matA12->getMat(x, y);
			val2 = matA22->getMat(x, y);
			sub1->setMat(x, y, val1-val2);
		}
	}

	for(x=0; x<sizeBRow2; x++) {
		for(y=0; y<sizeBColumn2; y++) {
			val3 = matB21->getMat(x, y);
			val4 = matB22->getMat(x, y);
			sum2->setMat(x, y, val3+val4);
		}
	}
	Strassen(sub1, sum2, matP7);

	//combine
	for(x=0; x<sizeARow2; x++) {
		for(y=0; y<sizeBColumn2; y++) {
			val1 = matP1->getMat(x, y);
			val2 = matP2->getMat(x, y);
			val3 = matP3->getMat(x, y);
			val4 = matP4->getMat(x, y);
			val5 = matP5->getMat(x, y);
			val6 = matP6->getMat(x, y);
			val7 = matP7->getMat(x, y);
			matC11->setMat(x, y, val1+val4-val5+val7);
			matC12->setMat(x, y, val3+val5);
			matC21->setMat(x, y, val2+val4);
			matC22->setMat(x, y, val1+val3-val2+val6);
		}
	}

	for(x=0; x<sizeARow1; x++) {
		for(y=0; y<sizeBColumn1; y++) {
			val1 = matC11->getMat(x, y);
			matC->setMat(x, y, val1);
		}
	}
	j = 0;

	for(x=0; x<sizeARow1; x++) {
		for(y=sizeBColumn1; y<matC->columns; y++) {
			val2 = matC12->getMat(x, j);
			matC->setMat(x, y, val2);
			j++;
		}
		j = 0;
	}
	i = 0;

	for(x=sizeARow1; x<matC->rows; x++) {
		for(y=0; y<sizeBColumn1; y++) {
			val3 = matC21->getMat(i, y);
			matC->setMat(x, y, val3);
		}
		i++;
	}
	i = 0; j = 0;

	for(x=sizeARow1; x<matC->rows; x++) {
		for(y = sizeBColumn1; y<matC->columns; y++) {
			val4 = matC22->getMat(i, j);
			matC->setMat(x, y, val4);
			j++;
		}
		i++;
		j = 0;
	}

	//release memory for temporal matrix
	matA11->freeMat(); matA12->freeMat(); matA21->freeMat(); matA22->freeMat();
	matB11->freeMat(); matB12->freeMat(); matB21->freeMat(); matB22->freeMat();
	matC11->freeMat(); matC12->freeMat(); matC21->freeMat(); matC22->freeMat();
	matP1->freeMat(); matP2->freeMat(); matP3->freeMat(); matP4->freeMat();
	matP5->freeMat(); matP6->freeMat(); matP7->freeMat();
	sum1->freeMat(); sum2->freeMat(); sub1->freeMat(); sub2->freeMat();
}

void Matrix::errorMat(const char *str) {
		perror(str);
		exit(EXIT_FAILURE);
}


int data1_[] = {
	12,	20,	16,	13,	5,	11,	10,	20,	8,	4,
	4,	11,	14,	10,	6,	2,	4,	18,	9,	15,
	2,	20,	1,	17,	14,	17,	19,	7,	10,	9,
	7,	6,	13,	2,	20,	5,	11,	18,	2,	11,
	19,	19,	8,	17,	12,	18,	10,	9,	3,	17,
	8,	3,	13,	10,	21,	17,	3,	18,	14,	3,
	6,	13,	4,	2,	15,	7,	2,	12,	5,	2,
	17,	15,	7,	5,	4,	6,	15,	6,	6,	20,
	2,	3,	4,	15,	10,	9,	2,	13,	8,	15,
	15,	15,	15,	16,	8,	4,	9,	7,	14,	3,
	17,	12,	8,	6,	13,	18,	6,	14,	13,	9,
	2,	20,	10,	8,	17,	14,	18,	7,	6,	2,
	10,	6,	18,	11,	2,	20,	16,	21,	12,	11,
	11,	18,	3,	17,	5,	9,	14,	11,	11,	6,
	4,	7,	21,	2,	21,	5,	16,	13,	4,	15,
	8,	14,	20,	19,	7,	19,	7,	14,	20,	2,
	20,	11,	9,	8,	11,	6,	2,	10,	15,	19,
	8,	8,	18,	20,	17,	6,	14,	9,	16,	17,
	12,	19,	15,	17,	9,	5,	3,	5,	3,	9,
	5,	6,	2,	14,	20,	9,	3,	20,	15,	5,
	2,	5,	5,	20,	5,	11,	10,	16,	21,	16,
	5,	14,	4,	19,	19,	11,	10,	5,	1,	7
 };

 int data2_[] = {
	7,	1,	8,	18,	5,	19,	5,	8,
	1,	17,	7,	7,	14,	14,	15,	14,
	14,	16,	17,	20,	4,	21,	5,	2,
	3,	3,	4,	6,	17,	21,	16,	2,
	21,	12,	9,	15,	2,	17,	4,	16,
	4,	8,	12,	15,	20,	17,	20,	3,
	9,	4,	14,	9,	12,	19,	12,	5,
	1,	5,	19,	6,	13,	15,	15,	14,
	12,	16,	19,	15,	13,	16,	4,	19,
	6,	21,	13,	18,	16,	10,	20,	17
 };

int main()
{
   Matrix *matA = new Matrix(22, 10, data1_);
   Matrix *matB = new Matrix(10, 8, data2_);
   Matrix *matC = new Matrix(22, 8);
   //matA->printMat();
   std::cout << std::endl;
   //matB->printMat();
   std::cout << std::endl;
   matC->Strassen(matA, matB, matC);
   matC->printMat();

   return 0;
}
