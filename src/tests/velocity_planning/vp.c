#include <stdio.h>
#include <math.h>
#include <raylib.h>
#include <util.h>
#include <vector.h>

#define WINDOW_WIDTH  1920
#define WINDOW_HEIGHT 1080

vec2_t vlim(double x, void* data) {
	return (vec2_t){
		0.25 * sin(10*x) + 0.3,
		2.5 * cos(10*x)
	};
}
double map(double x, double a, double b, double c, double d) {
	return c + (x - a) * (d - c) / (b - a);
}
void plot(vec2_t (*f)(double, void*), void* data) {
	int n = 200;
	double h = 1080;
	double w = 1920;
	double a = 0;
	double b = 1;
	for(int i = 0; i < n-1; i++) {
		double x0 = range(a, b, i+0, n);
		double x1 = range(a, b, i+1, n);
		vec2_t y0 = f(x0, data);
		vec2_t y1 = f(x1, data);
		
		if(y0.y * y1.y < 0) {
			DrawLine(
				map(x0, a, b, 0, w),
				map(y0.x, 0, 1, 0, h) - 10,
				map(x0, a, b, 0, w),
				map(y0.x, 0, 1, 0, h) + 10,
				RED
			);
		}

		DrawLine(
			map(x0, a, b, 0, w),
			map(y0.x, 0, 1, 0, h),
			map(x1, a, b, 0, w),
			map(y1.x, 0, 1, 0, h),
			WHITE
		);
	}
}

typedef struct alg_conf {
	vec2_t (*lim)(double, void*);
	double stride, t_max, v_max, a_max, j_max;
} alg_conf_t;

void alg(alg_conf conf) {
	//double (*lim)(double, void*) = conf.lim;	
}

int main(void) {
	InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Velocity Planning");	
	SetTargetFPS(144);

	alg_conf_t conf = {
		.lim    = vlim,
		.stride = 0.001,
		.t_max  = 1.0,
		.v_max  = 1.0,
		.a_max  = 0.1,
		.j_max  = 11.0
	};

	while(!WindowShouldClose()) {
		BeginDrawing();
		ClearBackground(BLACK);
		
		Vector2 mpos = GetMousePosition();
		int xpos = min(max(mpos.x, 0), WINDOW_WIDTH);
		int ypos = min(max(mpos.y, 0), WINDOW_HEIGHT);
		if(IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
			DrawLine(xpos, 0, xpos, WINDOW_HEIGHT, GRAY);
			conf.t_max = map((double)xpos, 0, WINDOW_WIDTH, 0, 1);	
			printf("t = %f\n", conf.t_max);
		} else if(IsKeyDown(KEY_ONE)) {
			DrawRectangle(WINDOW_WIDTH - 50, ypos, 50, WINDOW_HEIGHT - ypos, BLUE);
			conf.v_max = map((double)ypos, WINDOW_HEIGHT, 0, 0, 1);
			printf("v_max = %f\n", conf.v_max);
		} else if(IsKeyDown(KEY_TWO)) {
			DrawRectangle(WINDOW_WIDTH - 50, ypos, 50, WINDOW_HEIGHT - ypos, GREEN);		
			conf.a_max = map((double)ypos, WINDOW_HEIGHT, 0, 0, 10);
			printf("a_max = %f\n", conf.a_max);
		} else if(IsKeyDown(KEY_THREE)) {
			DrawRectangle(WINDOW_WIDTH - 50, ypos, 50, WINDOW_HEIGHT - ypos, RED);		
			conf.j_max = map((double)ypos, WINDOW_HEIGHT, 0, 0, 50);
			printf("j_max = %f\n", conf.j_max);
		}
		
		alg(conf);

		plot(vlim, NULL);

		EndDrawing();
	}
	
	CloseWindow();
	return 0;
}
