#define RL_CULL_DISTANCE_NEAR 0.01
#define RL_CULL_DISTANCE_FAR  1000.0

#include <stdio.h>
#include <math.h>
#include <raylib.h>
#include <raymath.h>
#include <rcamera.h>
#include <rlgl.h>
#include <imgui.h>
#include <rlImGui.h>

#include <opw.h>
#include <util.h>
#include <intrpl.h>
#include <eval.h>

#define WINDOW_WIDTH  1280
#define WINDOW_HEIGHT 720

#define COL_PURPLE (Color){127, 131, 255, 255}
#define COL_RED    (Color){255, 74,  152, 255}
#define COL_BLUE   (Color){10,  250, 250, 255}
#define COL_GRAY   (Color){127, 127, 127, 255}
#define COL_WHITE  (Color){255, 255, 255, 255}

#define CAMERA_MOVE_SPEED               0.09
#define CAMERA_ROTATION_SPEED           0.03
#define CAMERA_PAN_SPEED                0.2
#define CAMERA_MOUSE_MOVE_SENSITIVITY   0.003
#define CAMERA_MOUSE_SCROLL_SENSITIVITY 1.5
#define CAMERA_ORBITAL_SPEED            0.5

void update_camera(Camera *camera, int disable_keys, int disable_mouse) {
	Vector2 mouse_delta = GetMouseDelta();
	
	if(!disable_keys) {
		if(IsKeyDown(KEY_J)) { CameraPitch(camera, -CAMERA_ROTATION_SPEED, 0, 0, 0); }
		if(IsKeyDown(KEY_K)) { CameraPitch(camera, CAMERA_ROTATION_SPEED, 0, 0, 0); }
		if(IsKeyDown(KEY_L)) { CameraYaw(camera, -CAMERA_ROTATION_SPEED, 0); }
		if(IsKeyDown(KEY_H)) { CameraYaw(camera, CAMERA_ROTATION_SPEED, 0); }
		if(IsKeyDown(KEY_Q)) { CameraRoll(camera, -CAMERA_ROTATION_SPEED); }
		if(IsKeyDown(KEY_E)) { CameraRoll(camera, CAMERA_ROTATION_SPEED); }
		if(IsKeyDown(KEY_W)) { CameraMoveForward(camera, CAMERA_MOVE_SPEED, 0); }
		if(IsKeyDown(KEY_A)) { CameraMoveRight(camera, -CAMERA_MOVE_SPEED, 0); }
		if(IsKeyDown(KEY_S)) { CameraMoveForward(camera, -CAMERA_MOVE_SPEED, 0); }
		if(IsKeyDown(KEY_D)) { CameraMoveRight(camera, CAMERA_MOVE_SPEED, 0); }
		if(IsKeyDown(KEY_SPACE)) { CameraMoveUp(camera, CAMERA_MOVE_SPEED); }
		if(IsKeyDown(KEY_LEFT_CONTROL)) { CameraMoveUp(camera, -CAMERA_MOVE_SPEED); }
		if(IsKeyPressed(KEY_KP_SUBTRACT)) { CameraMoveToTarget(camera, 2.0f); }
		if(IsKeyPressed(KEY_KP_ADD))      { CameraMoveToTarget(camera, -2.0f); }
	}
	
	if(!disable_mouse) {
		if(IsKeyDown(KEY_LEFT_SHIFT) && IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
			if(mouse_delta.x > 0.0f) { CameraMoveRight(camera, CAMERA_PAN_SPEED, 0); }
			if(mouse_delta.x < 0.0f) { CameraMoveRight(camera, -CAMERA_PAN_SPEED, 0); }
			if(mouse_delta.y > 0.0f) { CameraMoveUp(camera, -CAMERA_PAN_SPEED); }
			if(mouse_delta.y < 0.0f) { CameraMoveUp(camera, CAMERA_PAN_SPEED); }
		} else if(IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
			CameraYaw(camera, -mouse_delta.x*CAMERA_MOUSE_MOVE_SENSITIVITY, 0);
			CameraPitch(camera, -mouse_delta.y*CAMERA_MOUSE_MOVE_SENSITIVITY, 0, 0, 0);
		}
		CameraMoveToTarget(camera, -GetMouseWheelMove());
	}
	
	return;
}

Camera3D setup_raylib(void) {
	SetConfigFlags(FLAG_WINDOW_RESIZABLE);
	SetConfigFlags(FLAG_MSAA_4X_HINT);
	InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "arm");

	Camera3D camera   = {0};
	camera.position   = (Vector3){10, 10, 10};
	camera.target     = (Vector3){0, 0, 0};
	camera.up         = (Vector3){0, 0, 1};
	camera.projection = CAMERA_PERSPECTIVE;
	camera.fovy       = 45.0;

	SetTargetFPS(60);

    rlImGuiBeginInitImGui();
    ImGui::StyleColorsDark();

	/*
    ImGuiIO& io = ImGui::GetIO();
    io.Fonts->AddFontFromFileTTF("font.ttf", 18);
	*/

    rlImGuiEndInitImGui();

	return camera;
}

Vector3 from_vec3t(vec3_t v) {
	return (Vector3){(float)v.x, (float)v.y, (float)v.z};
}

void draw_axes(void) {
	float l = 2;
	Vector3 origin = {0, 0, 0};
	Vector3 x_axis = {l, 0, 0};
	Vector3 y_axis = {0, l, 0};
	Vector3 z_axis = {0, 0, l};
	DrawLine3D(origin, x_axis, COL_RED);
	DrawLine3D(origin, y_axis, COL_PURPLE);
	DrawLine3D(origin, z_axis, COL_BLUE);
	rlDrawRenderBatchActive();
	return;
}

void raylib_matrix(double* mat, float* raymat) {
	for(int i = 0; i < 16; i++) {
		raymat[i] = mat[4*(i%4)+i/4];
	}
	return;
}

double degrees(double radians) {
	return radians*180/M_PI;
}

void draw_arm_joint(double angle) {
	double r = 0.3;
	DrawCylinderWires((Vector3){0, 0, 0}, 0.2, 0.2, 0, 20, COL_PURPLE);
	DrawCylinder((Vector3){0, -0.3, 0}, 0.0, 0.1, 0.3, 20, COL_PURPLE);
	DrawLine3D((Vector3){0, -0.3, 0}, (Vector3){0.3, -0.3, 0}, COL_BLUE);
	for(int i = 0; i < 20; i++) {
		double a1 = -(i/20.0)*angle;
		double a2 = -((i+1)/20.0)*angle;
		vec3_t p1 = (vec3_t){r*cos(a1), -0.3, r*sin(a1)};
		vec3_t p2 = (vec3_t){r*cos(a2), -0.3, r*sin(a2)};
		DrawLine3D(from_vec3t(p1), from_vec3t(p2), COL_BLUE);
	}
	rlPushMatrix();
	rlRotatef(degrees(angle), 0.0, 1.0, 0.0);
	DrawLine3D((Vector3){0, -0.3, 0}, (Vector3){0.3, -0.3, 0}, COL_BLUE);
	rlPopMatrix();
	return;
}

void draw_arm(kinstate_t* state, double* frames) {
	float t1[16], t2[16];

	rlPushMatrix();
	rlPushMatrix();
	rlRotatef(degrees(state->linkage[0].theta), 0.0, 0.0, 1.0);
	rlRotatef(90, 1.0, 0.0, 0.0);
	draw_arm_joint(-state->linkage[0].theta);
	rlPopMatrix();
	raylib_matrix(frames, t2);
	DrawLine3D((Vector3){0,0,0}, (Vector3){t2[12], t2[13], t2[14]}, COL_GRAY);
	for(int i = 0; i < 5; i++) {
		rlPushMatrix();
		raylib_matrix(frames+16*i, t1);
		rlMultMatrixf(t1);
		rlPushMatrix();
		rlRotatef(degrees(state->linkage[i+1].theta), 0.0, 0.0, 1.0);
		rlRotatef(90, 1.0, 0.0, 0.0);
		draw_arm_joint(-state->linkage[i+1].theta);
		rlPopMatrix();
		rlPopMatrix();
		raylib_matrix(frames+16*(i+1), t2);
		DrawLine3D((Vector3){t1[12], t1[13], t1[14]}, (Vector3){t2[12], t2[13], t2[14]}, COL_GRAY);
	}
	rlPopMatrix();	

	return;
}

vec3_t position_function1(ppoly_t ppoly, double t) {
	return eval_ppoly1(&ppoly, t).p;
}

vec3_t position_function3(ppoly_t ppoly, double t) {
	return eval_ppoly3(&ppoly, t).p;
}

vec3_t position_function5(ppoly_t ppoly, double t) {
	return eval_ppoly5(&ppoly, t).p;
}

vec3_t (*position_function)(ppoly_t, double) = &position_function3;

rotm_t orientation_function(double t) {
	return opw_euler(0, M_PI, 0);
}

void draw_path(ppoly_t ppoly) {
	for(int i = 1; i < 200; i++) {
		double t0 = range(0, 1, i-1, 200);
		double t1 = range(0, 1, i, 200);
		vec3_t p0 = position_function(ppoly, t0);
		vec3_t p1 = position_function(ppoly, t1);
		DrawLine3D(from_vec3t(p0), from_vec3t(p1), COL_GRAY);
	}
}

void draw_gui() {
	rlImGuiBegin();	
	//ImGui::ShowDemoWindow();
	rlImGuiEnd();
}

int main(void) {
	Camera3D camera = setup_raylib();	

	double sols[6*8];	
	double frames[6*4*4];
	int sol = 0;

	opw_t arm = {
		.c1 = 1, .c2 = 4, .c3 = 4, .c4 = 2,
		.a1 = 1, .a2 = -1,
		.b = 0
	};

	kinstate_t state = opw_kinstate(arm);

	vec3_t target_pos;
	rotm_t target_rot;

	ImGuiIO& io = ImGui::GetIO();

	double t = 0;

	/* testing */
	vec3_t points[5];
	points[0] = (vec3_t){4, 4, 1};
	points[1] = (vec3_t){-4, 4, -1};
	points[2] = (vec3_t){-4, -4, 1};
	points[3] = (vec3_t){4, -4, 0};
	points[4] = (vec3_t){5, -4, 1};
	
	vec3_t linear_mem0[EVAL_PPOLY1_MEMSIZE(5)];
	double linear_mem1[5];
	ppoly_t linear_ = eval_ppoly_init(linear_mem0, linear_mem1, 5);
	intrpl_ppoly1(points, 5, &linear_);

	vec3_t cubic_mem0[EVAL_PPOLY3_MEMSIZE(5)];
	double cubic_mem1[5];
	ppoly_t cubic_ = eval_ppoly_init(cubic_mem0, cubic_mem1, 5);
	int pivot3_mem[EVAL_PPOLY3_PIVOT_MEMSIZE(5)];
	double matrix3_mem[EVAL_PPOLY3_EQSYS_MEMSIZE(5)];
	double b3_mem[EVAL_PPOLY3_B_MEMSIZE(5)];

	vec3_t quintic_mem0[EVAL_PPOLY5_MEMSIZE(5)];
    double quintic_mem1[5];
    ppoly_t quintic_ = eval_ppoly_init(quintic_mem0, quintic_mem1, 5);
    int pivot5_mem[EVAL_PPOLY5_PIVOT_MEMSIZE(5)];
    double matrix5_mem[EVAL_PPOLY5_EQSYS_MEMSIZE(5)];
    double b5_mem[EVAL_PPOLY5_B_MEMSIZE(5)];

    luctx_t ctx3 = eval_ppoly_luctx_init(pivot3_mem, matrix3_mem, b3_mem);
    luctx_t ctx5 = eval_ppoly_luctx_init(pivot5_mem, matrix5_mem, b5_mem);

	intrpl_ppoly3(points, 5, NULL, NULL, &ctx3, &cubic_); 
	intrpl_ppoly5(points, 5, NULL, NULL, NULL, NULL, &ctx5, &quintic_); 

	ppoly_t active = cubic_;

	/* testing */
	
	int sign = 1;

	while(!WindowShouldClose()) {
		t += sign*0.002;
		if(t > 1.0) { t -= 0.002; sign *= -1; }
		if(t < 0.0) { t += 0.002; sign *= -1; }

		BeginDrawing();
		ClearBackground(BLACK);
		update_camera(&camera, io.WantCaptureKeyboard, io.WantCaptureMouse);
		BeginMode3D(camera);
		
		for(int i = 0; i < 5; i++) {
			DrawSphere(from_vec3t(points[i]), 0.1, COL_PURPLE);
		}
		
		draw_axes();

		KeyboardKey key = (KeyboardKey)GetKeyPressed();
		switch(key) {
			case KEY_ONE:   active = linear_; position_function = &position_function1; break;
			case KEY_TWO:   active = cubic_; position_function = &position_function3; break;
			case KEY_THREE: active = quintic_; position_function = &position_function5; break;
			default: break;
		}

		/* main logic */	
		target_pos = position_function(active, t);
		target_rot = orientation_function(t);	
		draw_path(active);
		opw_inverse(arm, target_pos, target_rot, sols);	
		opw_forward(state, frames);
		opw_update_kinstate(&state, sols + 6*sol);
		draw_arm(&state, frames);	
		/* main logic */

		EndMode3D();
		
		draw_gui();

		EndDrawing();
	}

	CloseWindow();

	return 0;
}
