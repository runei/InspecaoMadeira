#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> /* define inteiros de tamanho específico */
#include <string.h>
#include <dirent.h>
#include <float.h>
#include <math.h>
#include <gperftools/profiler.h>
#include <omp.h>

#define K_WOOD 128
#define LEFT -1
#define RIGHT 1
#define UP 1
#define DOWN -1
#define WHITE 0xFFFFFF
#define BLACK 0x000000
#define BLUE 0x0000FF
#define GREEN 0x00FF00
#define BLACK_BACKGOURND_WOOD 0x0C0C0C
#define NTHREADS 2
#define H_SIZE 512
#define W_SIZE 488

#pragma pack(push, 1)
typedef struct {
	uint16_t type;
	uint32_t size;
	uint16_t reserved1;
	uint16_t reserved2;
	uint32_t offset;
	uint32_t header_size;
	int32_t width;
	int32_t height;
	uint16_t planes;
	uint16_t bits;
	uint32_t compression;
	uint32_t imagesize;
	int32_t xresolution;
	int32_t yresolution;
	uint32_t ncolours;
	uint32_t importantcolours;
} BitmapHeader;
#pragma pack(pop) /* restaura comportamento do compilador */

typedef struct {
	BitmapHeader bmp_header;
	unsigned int matrix[H_SIZE][W_SIZE];
} Bitmap;

typedef struct {
	int row;
	int col;
} Point;

typedef struct {
	int size;
	Point *nodes;
} Nodes;

typedef struct {
	int stack_pointer;
	Point *stack;
} Stack;

uint8_t getBlue(unsigned int color)
{
	return 0xff & color;
}

uint8_t getGreen(unsigned int color)
{
	return (color >> 8) & 0xff;
}

uint8_t getRed(unsigned int color)
{
	return (color >> 16) & 0xff;
}

uint8_t getGray(unsigned int color)
{
	return (getRed(color) + getGreen(color) + getBlue(color)) / 3;
}

unsigned int setColor(int red, int green, int blue)
{
	return (red << 16) | (green << 8) | blue;
}

double getValue(unsigned int color)
{
	return getGray(color) / 255.0 * 100.0;
}

void getIniEndParallelFor(const int old_init, const int old_end, int *new_init, int *new_end)
{
	const int nthreads = omp_get_num_threads();
	const int part = (old_end - old_init) / nthreads;
	const int id = omp_get_thread_num();
	*new_init = (part * id) + old_init;

	if (id == nthreads - 1)
	{
		*new_end = old_end;
	}
	else
	{
		*new_end = *new_init + part;
	}
}

void push(Stack *stack, Point val)
{
	if (stack == NULL)
	{
		perror("Stack = NULL");
		return;
	}
	++stack->stack_pointer;
	Point *p = (Point *)realloc(stack->stack, (stack->stack_pointer+1) * sizeof(Point));
	if (p != NULL)
	{
		stack->stack = p;
	}
	stack->stack[stack->stack_pointer] = val;
}

Point pop(Stack *stack)
{
	Point p = stack->stack[stack->stack_pointer];
	--stack->stack_pointer;
	return p;
}

int max(int a, int b)
{
	if (a > b)
	{
		return a;
	}
	return b;
}

int min(int a, int b)
{
	if (a < b)
	{
		return a;
	}
	return b;
}

unsigned int **sobelFiltering(Bitmap *bmp)
{
	unsigned int **result = (unsigned int **) malloc(sizeof(unsigned int *) * bmp->bmp_header.height);
	for (int i = 0; i < bmp->bmp_header.height; i++) {
		result[i] = (unsigned int *) malloc(sizeof(unsigned int) * bmp->bmp_header.width);
	}
  /* Definition of Sobel filter in horizontal direction */
	int weight_x[3][3] = {{ -1,  0,  1 },
						{ -2,  0,  2 },
						{ -1,  0,  1 }};
	int weight_y[3][3] = {{ 1,  2,  1 },
						{ 0,  0,  0 },
						{ -1,  -2,  -1 }};

  /* Maximum values calculation after filtering*/
	for (int i = 1; i < bmp->bmp_header.height - 1; ++i)
	{
		for (int j = 1; j < bmp->bmp_header.width - 1; ++j)
		{
			int gx = 0;
			int gy = 0;
			for (int k = 0; k < 3; ++k)
			{
				for (int l = 0; l < 3; ++l)
				{
					int val = getGray(bmp->matrix[i+k-1][j+l-1]);
					gx += val * weight_x[k][l];
					gy += val * weight_y[k][l];
				}
			}

			// printf("%u\n", val);
			int val = (int) sqrt(pow(gx, 2) + pow(gy, 2));
			if (val < 220)
			{
				val = 255;
			}
			else //if (val < 0)
			{
				val = 0;
			}
			result[i][j] = setColor(val, val, val);
			// result[i][j] = abs(gx) + abs(gy);
			// printf("%u\n", setColor(K_WOOD, K_WOOD, K_WOOD));
			/*if (result[i][j] < setColor(K_WOOD, K_WOOD, K_WOOD))
			{
				result[i][j] = BLACK;
			}
			else
			{
				result[i][j] = WHITE;
			}*/
			// result[i][j] = sqrt(pow(gx, 2) + pow(gy, 2));
// exit(0);
		}
	}
	return result;
}

void loadBitmap(const char *file_name, Bitmap *bitmap)
{
	FILE *file;

	file = fopen(file_name, "rb");

	if (file == NULL)
	{
		printf("Erro ao abrir a imagem %s\n", file_name);
		perror("");
		exit(0);
	}

	if (!fread(&bitmap->bmp_header, sizeof(BitmapHeader), 1, file))
	{
		perror("Erro ao ler bitmap header");
		exit(0);
	}

	//ler o arquivo
	const long fpoint = ftell(file);
	fseek(file, 0, SEEK_END);
	const long fsize = ftell(file) - fpoint;
	fseek(file, fpoint, SEEK_SET);

	char *image = malloc(fsize + 1);
	if (!fread(image, fsize, 1, file))
	{
		perror("Erro ao ler bitmap");
		exit(0);
	}
	fclose(file);

	int index = 0;

	/* leitura */
	for (int i = 0; i < bitmap->bmp_header.height; i++) {
		for (int j = 0; j < bitmap->bmp_header.width; j++) {
			const unsigned int color = setColor(image[index + 2], image[index + 1], image[index]);
			const unsigned int val = getGray(color);
			bitmap->matrix[i][j] = setColor(val, val, val);
			index += 3;
		}
	}

	free(image);
	// return matrix;
}

void saveBitmap(char *fileName, BitmapHeader bH, unsigned int matrix[H_SIZE][W_SIZE])
{
	FILE *nova = fopen(fileName, "wb");

	fwrite(&bH, sizeof(BitmapHeader), 1, nova);

	/* escrita */
	for (int i = 0; i < bH.height; i++) {
		for (int j = 0; j < bH.width; j++) {
			fwrite(&matrix[i][j], bH.bits / 8, 1, nova); // 1byte = 8bits
		}
	}

	fclose(nova);
}

void printBitmap(const BitmapHeader bH, const unsigned int matrix[H_SIZE][W_SIZE])
{
	for (int i = 0; i < bH.height; i++) {
		for (int j = 0; j < bH.width; j++) {
			printf("%d ", matrix[i][j]);
		}
		printf("\n");
	}
}

int invertPosMatrix(const int height, const int val)
{
	return height - val - 1;
}

int getLinhaInicial(const BitmapHeader bH, const unsigned int matrix[H_SIZE][W_SIZE])
{
	//nao considera imagens vazias
	int i = bH.height - 1;
	while (getValue(matrix[i][10]) < 10.0)
	{
		--i;
	}
	return invertPosMatrix(bH.height, i);
}

int getLinhaFinal(const BitmapHeader bH, const unsigned int matrix[H_SIZE][W_SIZE])
{
	//nao considera imagens vazias
	int i = 0;
	while (getValue(matrix[i][10]) < 10.0)
	{
		++i;
	}
	return invertPosMatrix(bH.height, i);
}

int checkErosion(int i, int j, unsigned int matrix[H_SIZE][W_SIZE])
{
	for (int k = -1; k < 2; ++k)
	{
		for (int l = -1; l < 2; ++l)
		{
			if (matrix[i-k][j-l] == WHITE)
			{
				return 1;
			}
		}
	}
	return 0;
}

int checkDilation(const int i, const int j, const int size, unsigned int matrix[H_SIZE][W_SIZE])
{
	const int ini = -size;
	const int fim = size + 1;

	int result = 0;

	for (int k = ini; k < fim; ++k)
	{
		for (int l = ini; l < fim; ++l)
		{
			result |= matrix[i+k][j+l] == BLACK;
			// printf("matrix[i+%d][j+%d] == BLACK || ", k, l);
		}
	}
	return result;
}

void dilation(const int init_line, const int end_line, Bitmap *bitmap)
{
	unsigned int flag_matrix[H_SIZE][W_SIZE];
	char if_matrix[H_SIZE][W_SIZE];
	memset(if_matrix, 0, H_SIZE * W_SIZE);
	const int size = end_line - init_line;
	// flag_matrix = (unsigned int **)malloc(sizeof(unsigned int *) * size);
	//matriz bitmap é invertida
	const int temp_init = invertPosMatrix(bitmap->bmp_header.height, end_line);
	const int temp_end = invertPosMatrix(bitmap->bmp_header.height, init_line);

	static const int size_dilation = 2;

	int id, i, j, part, ini, fim, nthreads;

	const int max_height = temp_end - size_dilation;
	const int max_width = bitmap->bmp_header.width - size_dilation;

/*	for (i = temp_init + size_dilation; i < max_height; ++i)
	{
		for (j = size_dilation; j < max_width; ++j)
		{
			const char test = bitmap->matrix[i][j] == BLACK;
			if_matrix[i-2][j] |= test;
			if_matrix[i-1][j] |= test;
			if_matrix[i][j-1] |= test;
			if_matrix[i][j-2] |= test;
			if_matrix[i][j] |= test;
		}
		if_matrix[i-2][max_width-2] |= bitmap->matrix[i-2][max_width] == BLACK;
		if_matrix[i-2][max_width-2] |= bitmap->matrix[i-2][max_width-1] == BLACK;
		if_matrix[i-1][max_width-2] |= bitmap->matrix[i-1][max_width] == BLACK;
		if_matrix[i-1][max_width-2] |= bitmap->matrix[i-1][max_width-1] == BLACK;
	}

	for (j = size_dilation; j < max_width; ++j)
	{
		if_matrix[max_height-2][j] |= bitmap->matrix[max_height][j] == BLACK;
		if_matrix[max_height-2][j] |= bitmap->matrix[max_height-1][j] == BLACK;
	}
*/	// if_matrix[i][j-1] |= test;
	// if_matrix[i][j-2] |= test;

	#pragma omp parallel private(i, j, ini, fim)
	{
// printf("%d\t%d\n", omp_get_num_threads(), omp_get_thread_num());
		getIniEndParallelFor(temp_init, temp_end, &ini, &fim);

		for (i = ini; i < fim; ++i)
		// for (i = temp_init; i < temp_end; ++i)
		{
			for (j = 0; j < bitmap->bmp_header.width; ++j)
			{
				if (i < temp_init + size_dilation || i >= max_height || j < size_dilation || j >= max_width)
				{
					flag_matrix[i][j] = BLACK;
				}
				// else if (!if_matrix[i][j])
				// else if (!checkDilation(i, j, size_dilation, bitmap->matrix))
				else if (!(bitmap->matrix[i-2][j-2] == BLACK || bitmap->matrix[i-2][j-1] == BLACK || bitmap->matrix[i-2][j] == BLACK || bitmap->matrix[i-2][j+1] == BLACK || bitmap->matrix[i-2][j+2] == BLACK || bitmap->matrix[i-1][j-2] == BLACK || bitmap->matrix[i-1][j-1] == BLACK || bitmap->matrix[i-1][j] == BLACK || bitmap->matrix[i-1][j+1] == BLACK || bitmap->matrix[i-1][j+2] == BLACK || bitmap->matrix[i][j-2] == BLACK || bitmap->matrix[i][j-1] == BLACK || bitmap->matrix[i][j] == BLACK || bitmap->matrix[i][j+1] == BLACK || bitmap->matrix[i][j+2] == BLACK || bitmap->matrix[i+1][j-2] == BLACK || bitmap->matrix[i+1][j-1] == BLACK || bitmap->matrix[i+1][j] == BLACK || bitmap->matrix[i+1][j+1] == BLACK || bitmap->matrix[i+1][j+2] == BLACK || bitmap->matrix[i+2][j-2] == BLACK || bitmap->matrix[i+2][j-1] == BLACK || bitmap->matrix[i+2][j] == BLACK || bitmap->matrix[i+2][j+1] == BLACK || bitmap->matrix[i+2][j+2] == BLACK))
				{
					// exit(0);
					flag_matrix[i][j] = WHITE;
				}
				else
				{
					flag_matrix[i][j] = BLACK;
				}
			}
		}
	// }
	// }
		#pragma omp barrier
		// #pragma omp taskwait

		// for (i = temp_init; i < temp_end; ++i)
		for (i = ini; i < fim; ++i)
		{
			for (j = 0; j < bitmap->bmp_header.width; ++j)
			{
				bitmap->matrix[i][j] = flag_matrix[i][j];
			}
		}
	}
}

void erosion(const int init_line, const int end_line, Bitmap *bitmap)
{
	unsigned int flag_matrix[H_SIZE][W_SIZE];
	//matriz bitmap é invertida
	const int temp_init = invertPosMatrix(bitmap->bmp_header.height, end_line);
	const int temp_end = invertPosMatrix(bitmap->bmp_header.height, init_line);

	for (int i = temp_init; i < temp_end; ++i)
	{
		for (int j = 0; j < bitmap->bmp_header.width; ++j)
		{
			if (i == temp_init || i == temp_end - 1 || j == 0 || j == bitmap->bmp_header.width - 1 || checkErosion(i, j, bitmap->matrix))
			{
				flag_matrix[i - temp_init][j] = WHITE;
			}
			else
			{
				flag_matrix[i - temp_init][j] = BLACK;
			}
		}
	}

	for (int i = temp_init; i < temp_end; ++i)
	{
		for (int j = 0; j < bitmap->bmp_header.width; ++j)
		{
			bitmap->matrix[i][j] = flag_matrix[i - temp_init][j];
		}
	}

}

int fill(int temp_init, Point p, Bitmap *bmp, char flag_matrix[H_SIZE][W_SIZE], char side)
{
	int result = p.col;
	do
	{
		flag_matrix[p.row-temp_init][result] = 1;
		// bmp->matrix[p.row][result] = color;
		result += side;
	} while (result < bmp->bmp_header.width && result >= 0 && bmp->matrix[p.row][result] == BLACK);
	return result - side;
}

void findStartPoints(int xleft, int xright, int col, unsigned int *row, Stack *stack, char *row_flag_matrix)
{
	char new_point = 1;
	int i = xleft;
	while (i <= xright)
	{
		if (new_point && !row_flag_matrix[i] && row[i] == BLACK)
		{
			Point p = {col, i};
			push(stack, p);
			new_point = 0;
		}
		else if (!new_point && row[i] == WHITE)
		{
			new_point = 1;
		}
		++i;
	}
}

Nodes getWoodNodes(const int init_line, const int end_line, Bitmap *bitmap)
{
	char flag_matrix[H_SIZE][W_SIZE];
	int size = end_line - init_line;

	//matriz bitmap é invertida
	const int temp_init = invertPosMatrix(bitmap->bmp_header.height, end_line);
	const int temp_end = invertPosMatrix(bitmap->bmp_header.height, init_line);

	Nodes nodes = {0, NULL};

	memset(flag_matrix, 0, H_SIZE * W_SIZE);

	for (int i = temp_init; i < temp_end; ++i) // nao conta contornos das linhas
	{
		for (int j = 0; j < bitmap->bmp_header.width; ++j)
		{
			if (!flag_matrix[i-temp_init][j])
			{
				flag_matrix[i-temp_init][j] = 1;
				if (bitmap->matrix[i][j] == BLACK)
				{

					int max_col, min_col, max_row, min_row;
					max_col = min_col = j;
					max_row = min_row = i;
					Stack stack = {0, NULL};

					Point p = {i, j};
					push(&stack, p);
					while (stack.stack_pointer)
					{
						p = pop(&stack);
						max_col = max(max_col, p.col);
						min_col = min(min_col, p.col);
						max_row = max(max_row, p.row);
						min_row = min(min_row, p.row);

						int xleft, xright;

						xleft = fill(temp_init, p, bitmap, flag_matrix, LEFT);
						xright = fill(temp_init, p, bitmap, flag_matrix, RIGHT);
						if (p.row-temp_init > 0)
						{
							findStartPoints(xleft, xright, p.row-1, bitmap->matrix[p.row-1], &stack, flag_matrix[p.row-temp_init-1]);
						}
						if (p.row-temp_init < size - 1)
						{
							findStartPoints(xleft, xright, p.row+1, bitmap->matrix[p.row+1], &stack, flag_matrix[p.row-temp_init+1]);
						}
					}

					if (min_row > temp_init && max_row < temp_end - 1)
					{
						++nodes.size;
						nodes.nodes = (Point *)realloc(nodes.nodes, nodes.size * sizeof(Point));
						nodes.nodes[nodes.size - 1].col = (max_col + min_col) / 2;
						nodes.nodes[nodes.size - 1].row = invertPosMatrix(bitmap->bmp_header.height, (max_row + min_row) / 2);
					}
					free(stack.stack);
/*
					char name[50];
					snprintf(name,50,"testes/nova_imagem%d_%d.bmp", i, j);
					saveBitmap(name, bitmap.bmp_header, bitmap.matrix);*/
				}
			}
		}
	}

	return nodes;
}

void contrast(const int init_line, const int end_line, Bitmap *bitmap)
{
	int i, j, ini, end, value, new_val;
	int max_val = 255, min_val = 0;

	const int temp_init = invertPosMatrix(bitmap->bmp_header.height, end_line);
	const int temp_end = invertPosMatrix(bitmap->bmp_header.height, init_line);

	int histogram[256];
	memset(histogram, 0, sizeof(int) * 256);

		// #pragma omp critical(histogram)
		for (i = temp_init; i < temp_end; i++) {
		// for (i = ini; i < end; i++) {
			for (int j = 0; j < bitmap->bmp_header.width; j++) {
				const int color = getGray(bitmap->matrix[i][j]);
				// #pragma omp critical
				{
					histogram[color]++;
				}
			}
		}

		// #pragma omp barrier

		// #pragma omp single nowait
		{
			i = min_val;
			while (histogram[i] == 0)
			{
				min_val++;
				i++;
			}
		}

		// #pragma omp single
		{
			i = max_val;
			while (histogram[i] == 0)
			{
				max_val--;
				i--;
			}
		}

	#pragma omp parallel private(i, j, ini, end) shared(max_val, min_val)
	{

		getIniEndParallelFor(temp_init, temp_end, &ini, &end);


		for (i = ini; i < end; ++i) {
		// for (int i = temp_init; i < temp_end; i++) {
			for (j = 0; j < bitmap->bmp_header.width; j++) {
				value = getGray(bitmap->matrix[i][j]);

				new_val = (255 * (value - min_val) / (max_val - min_val));

				if (new_val > 255)
				{
					new_val = 255;
				}
				else if (new_val < 0)
				{
					new_val = 0;
				}
				bitmap->matrix[i][j] = setColor(new_val, new_val, new_val);
			}
		}
	}
}

int otsuThresholder(Bitmap *bitmap)
{
	int histogram[256];
	for(int i = 0; i < 256; ++i)
	{
		histogram[i] = 0;
	}

	for (int i = 0; i < bitmap->bmp_header.height; i++) {
		for (int j = 0; j < bitmap->bmp_header.width; j++) {
			int color = getGray(bitmap->matrix[i][j]);
			histogram[color]++;
		}
	}

/*	FILE *file = fopen("teste.txt", "w+");
	for(int i = 0; i < 256; ++i)
	{
		int x = histogram[i];
		// if (x > 20) x = 0;
		fprintf(file, "%d %d\n", i, x);
	}
	fclose(file);
*/
	int total = bitmap->bmp_header.height * bitmap->bmp_header.width;
	// int total = (temp_end - temp_init) * bitmap->bmp_header.width;

	double sum = 0;
	for (int i = 0; i < 256 ; i++)
	{
		sum += i * histogram[i];
	}

	double sum_b = 0;
	int wb = 0;
	int wf = 0;

	double var_max = 0;
	int threshold = 0;

	for (int  i = 0; i < 256; ++i)
	{
		wb += histogram[i];
		if (wb == 0) continue;

		wf = total - wb;
		if (wf == 0) break;

		sum_b += (double) (i * histogram[i]);

		double mb = sum_b / wb;
		double mf = (sum - sum_b) / wf;

		double var_between = (double) wb * (double) wf * (mb - mf) * (mb - mf);

		if (var_between > var_max) {
			var_max = var_between;
			threshold = i;
		}
	}
	return threshold;
}

void binarizeBitmap(int init_line, int end_line, Bitmap *bitmap, int threshold)
{
	int temp_init = invertPosMatrix(bitmap->bmp_header.height, end_line);
	int temp_end = invertPosMatrix(bitmap->bmp_header.height, init_line);

	for (int i = temp_init; i < temp_end; i++) {
		for (int j = 0; j < bitmap->bmp_header.width; j++) {
			if (getGray(bitmap->matrix[i][j]) > threshold)
			{
				bitmap->matrix[i][j] = WHITE;
			}
			else
			{
				bitmap->matrix[i][j] = BLACK;
			}
		}
	}
}

double test_program()
{
	FILE *result = fopen("result.txt", "r");

	int total_files = 0;
	int correct_files = 0;

	char *buffer_result = NULL;
	size_t buffsize_result = 0;
	ssize_t nread_result;
	while ((nread_result = getline(&buffer_result, &buffsize_result, result)) != -1)
	{
		char img_name[15];
		int pos_init, pos_end, n_nodes_result;
		sscanf(buffer_result, "%s %d %d %d", img_name, &pos_init, &pos_end, &n_nodes_result);

		FILE *test = fopen("test.txt", "r");

		char *buffer_test = NULL;
		size_t buffsize_test = 0;
		ssize_t nread_test;

		while ((nread_test = getline(&buffer_test, &buffsize_test, test)) != -1)
		{
			char img_num[10];
			int n_nodes_test;

			sscanf(buffer_test, "%s %d", img_num, &n_nodes_test);

			if (strstr(img_name, img_num) != NULL)
			{
				++total_files;
				if (n_nodes_result == n_nodes_test)
				{
					++correct_files;
				}
				break;
			}
		}

		fclose(test);
		free(buffer_test);
	}

	free(buffer_result);
	fclose(result);

	// printf("%d%d\n", total_files, correct_files);
	if (total_files > 0)
	{
		return (double) correct_files / (double) total_files * 100.0;
	}
	else
	{
		return 0.0;
	}
}

int main(int argc, char *argv[])
{
	// ProfilerStart("prof.log");
	// test();return 0;

	if (argc != 2)
	{
		perror("erro");
		return 0;
	}

	char *dir_name = argv[1];

	omp_set_num_threads(NTHREADS);
	omp_set_nested(1);
	// omp_set_dynamic(1);

	// #pragma omp parallel
	{
		// #pragma omp single
		{

	for (unsigned int i = 28; i < 29; ++i)
	// for (unsigned int i = 49; i < 50; ++i)
	{

		DIR *dir = opendir(dir_name);
		if (dir == NULL)
		{
			perror("Diretório 'validos' não encontrado");
			exit(0);
		}

		FILE *output = fopen("result.txt", "w+");

		struct dirent *ent;
		while ((ent = readdir(dir)) != NULL)
		{
			char *ext = strrchr(ent->d_name, '.');
			// if (strcmp(ent->d_name, "st1283.bmp") == 0)
			// if (strcmp(ent->d_name, "st1012.bmp") == 0)
			// if (strcmp(ent->d_name, "st1355.bmp") == 0)
			// if (strcmp(ent->d_name, "st1049.bmp") == 0)
			// if (strcmp(ent->d_name, "st1551.bmp") == 0)
			if (ext != NULL && strcmp(ext, ".bmp") == 0)
			{
				Bitmap bitmap;
				char file_name[300];
				snprintf(file_name, sizeof(file_name), "%s%s", dir_name, ent->d_name);
				loadBitmap(file_name, &bitmap);

				// saveBitmap("nova_imagem.bmp", bitmap.bmp_header, bitmap.matrix);
				// exit(0);

				int init_line = getLinhaInicial(bitmap.bmp_header, bitmap.matrix);
				int end_line = getLinhaFinal(bitmap.bmp_header, bitmap.matrix);

				contrast(init_line, end_line, &bitmap);

				int threshold = otsuThresholder(&bitmap);

				binarizeBitmap(init_line, end_line, &bitmap, threshold);
				unsigned char test = i;
				while (test)
				{
					if (test & 1)
					{
						dilation(init_line, end_line, &bitmap);
					}
					else
					{
						erosion(init_line, end_line, &bitmap);
					}
					test = test >> 1;
				}

				// saveBitmap("nova_imagem2.bmp", bitmap.bmp_header, bitmap.matrix);

// exit(0);

				Nodes nodes;
				do
				{
					nodes = getWoodNodes(init_line, end_line, &bitmap);
					// saveBitmap(ent->d_name, bitmap.bmp_header, bitmap.matrix);
					dilation(init_line, end_line, &bitmap);
				} while (nodes.size > 3);


				/*printf("%s %d %d %d", ent->d_name, init_line, end_line, nodes.size);
				for (int i = 0; i < nodes.size; ++i)
				{
					printf(" %d %d", nodes.nodes[i].row, nodes.nodes[i].col);
				}
				printf("\n");*/
				fprintf(output, "%s %d %d %d", ent->d_name, init_line, end_line, nodes.size);
				for (int i = 0; i < nodes.size; ++i)
				{
					fprintf(output, " %d %d", nodes.nodes[i].row, nodes.nodes[i].col);
				}
				fprintf(output, "\n");
				free(nodes.nodes);
			}
		}
		fclose(output);
		closedir(dir);
//		printf("Finished program. Init Test.\n");
		double perc = test_program();
		printf("%f\t%d\n", perc, i);
	}
	}}
	// ProfilerStop();

	return 0;
}
