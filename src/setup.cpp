#include "setup.hpp"

void main_setup() { // 带有阵风循环的城市风环境模拟; required extensions in defines.hpp: D3Q27, TRT, EQUILIBRIUM_BOUNDARIES, SUBGRID, GRAPHICS
	// ################################################################## define simulation box size, viscosity and volume force ###################################################################
	const uint3 lbm_N = resolution(float3(9.0f, 12.0f, 2.0f), 8000u); //定义模拟盒长宽高的比例为Rhino中的BoundingBox；显存占用约 8GB
	const float Re = 50000.0f;
	const float u_base = 0.08f;  // 基础风速 (lu)
	const float u_gust = 0.12f;  // 阵风风速 (lu)
	const ulong lbm_T = 5000ull; // 总模拟时间步数
	
	// ========== 物理单位转换 ==========
	const float U_ref_real = 5.0f;  // 赫尔辛基冬季参考风速 (m/s)
	const float S_u = U_ref_real / u_base;  // 速度缩放因子 = 62.5
	
	print_info("=== VELOCITY SCALING ===");
	print_info("Reference wind speed: " + to_string(U_ref_real, 1u) + " m/s");
	print_info("Base velocity (lu): " + to_string(u_base));
	print_info("Gust velocity (lu): " + to_string(u_gust));
	print_info("Scaling factor S_u: " + to_string(S_u, 1u));
	print_info("Physical base wind: " + to_string(u_base * S_u, 1u) + " m/s");
	print_info("Physical gust wind: " + to_string(u_gust * S_u, 1u) + " m/s");
	
	LBM lbm(lbm_N, units.nu_from_Re(Re, (float)lbm_N.z * 0.3f, u_base));

	// ###################################################################################### 定义几何体 ######################################################################################
	Mesh* building = read_stl(get_exe_path() + "../stl/test_terrace_2.stl");

	print_info("=== ORIGINAL MESH INFO ===");
	print_info("Size (W x D x H): " + to_string(building->get_bounding_box_size().x) + " x " +
		to_string(building->get_bounding_box_size().y) + " x " +
		to_string(building->get_bounding_box_size().z));

	const float target_height = 0.8f * (float)lbm_N.z; // 目标建筑高度为模拟盒Z方向的80%
	const float scale = target_height / building->get_bounding_box_size().z;
	building->scale(scale);

	// ========== 旋转建筑物（模拟斜风向） ==========
	// 绕 Z 轴旋转，正值为逆时针，负值为顺时针
	// 旋转 -15° 相当于风从 Y 轴方向偏转 -15° 吹向建筑
	const float rotation_angle = 12.9f; // 旋转角度（度）
	const float3x3 rotation = float3x3(float3(0, 0, 1), radians(rotation_angle));
	building->rotate(rotation);

	const float3 target_position = float3(
		0.5f * (float)lbm_N.x,
		0.4f * (float)lbm_N.y,
		1.0f + building->get_bounding_box_size().z * 0.5f
	);
	building->translate(target_position - building->get_bounding_box_center());

	print_info("=== AFTER SCALING ===");
	print_info("Scale factor: " + to_string(scale));

	print_info("=== VOXELIZING MESH ===");
	lbm.voxelize_mesh_on_device(building);

	lbm.flags.read_from_device();
	uint building_count = 0u;
	const uint Nx = lbm.get_Nx(), Ny = lbm.get_Ny(), Nz = lbm.get_Nz();
	for (ulong n = 0ull; n < lbm.get_N(); n++) {
		if (lbm.flags[n] == TYPE_S) {
			uint x = 0u, y = 0u, z = 0u;
			lbm.coordinates(n, x, y, z);
			if (z > 0u) building_count++;
		}
	}
	print_info("Building cells: " + to_string(building_count));

	// ###################################################################################### 设置边界条件和初始速度 ######################################################################################
		parallel_for(lbm.get_N(), [&](ulong n) {
		uint x = 0u, y = 0u, z = 0u;
		lbm.coordinates(n, x, y, z);

		if (z == 0u && lbm.flags[n] != TYPE_S) {
			lbm.flags[n] = TYPE_S;
		}
		if (lbm.flags[n] != TYPE_S) {
			lbm.u.x[n] = u_base;  // ← 改为 X 方向
		}
		if (x == 0u && lbm.flags[n] != TYPE_S) {  // ← 入口改为 x=0
			lbm.flags[n] = TYPE_E;
			lbm.u.x[n] = u_base;
		}
		if (x == Nx - 1u && lbm.flags[n] != TYPE_S) {  // ← 出口改为 x=Nx-1
			lbm.flags[n] = TYPE_E;
		}
		if ((y == 0u || y == Ny - 1u || z == Nz - 1u) && lbm.flags[n] != TYPE_S) {  // ← 侧边界
			lbm.flags[n] = TYPE_E;
		}
		});

	lbm.flags.write_to_device();
	lbm.u.write_to_device();

	// ###################################################################################### 运行模拟并渲染 ######################################################################################
	// ========== 可视化模式选择（城市风模拟推荐） ==========

	// 方案1: 固体表面 + Q准则涡旋（推荐，显示湍流结构）
	//lbm.graphics.visualization_modes = VIS_FLAG_SURFACE | VIS_Q_CRITERION;

	// 方案2: 只显示固体表面（最干净）
	lbm.graphics.visualization_modes = VIS_FLAG_SURFACE;

	// 方案3: 固体表面 + 流线
	// lbm.graphics.visualization_modes = VIS_FLAG_SURFACE | VIS_STREAMLINES;

	// 方案4: 固体表面 + 速度场切片
	//lbm.graphics.visualization_modes = VIS_FLAG_SURFACE | VIS_STREAMLINES;
	//lbm.graphics.slice_mode = 1; // x方向切片（垂面）
	//lbm.graphics.slice_x = 346; // 切片位置
	//lbm.graphics.field_mode = 0;   // 0=速度, 1=密度, 2=温度

	print_info("Grid size: Nx=" + to_string(Nx) + ", Ny=" + to_string(Ny) + ", Nz=" + to_string(Nz));

	print_info("=== STARTING GUST SIMULATION ===");

	lbm.run(0u, lbm_T); // 初始化模拟

	uint cycle = 0u;
	while (lbm.get_t() < lbm_T) {
		print_info("=== CYCLE " + to_string(cycle) + ": Normal wind ===");

		// ========== 正常风速 ==========
		lbm.u.read_from_device();
		for (ulong n = 0ull; n < lbm.get_N(); n++) {
			uint x = 0u, y = 0u, z = 0u;
			lbm.coordinates(n, x, y, z);
			if (x == 0u && lbm.flags[n] == TYPE_E) {
				lbm.u.x[n] = u_base;
			}
		}
		lbm.u.write_to_device();
		lbm.run(1000u, lbm_T);

		// ========== 速度统计（Normal wind）—— 仅监控视野内区域 ==========
		lbm.u.read_from_device();
		float u_max_normal = 0.0f;
		bool danger_normal = false;

		// 监控区域定义（使用 defines.hpp 中的宏）
		const uint x_min = (uint)(MONITOR_X_MIN * Nx), x_max = (uint)(MONITOR_X_MAX * Nx);
		const uint y_min = (uint)(MONITOR_Y_MIN * Ny), y_max = (uint)(MONITOR_Y_MAX * Ny);
		const uint z_min = MONITOR_Z_MIN_CELLS, z_max = (uint)(MONITOR_Z_MAX * Nz);

		for (ulong n = 0ull; n < lbm.get_N(); n++) {
			if (lbm.flags[n] != TYPE_S) {
				uint x = 0u, y = 0u, z = 0u;
				lbm.coordinates(n, x, y, z);
				const float u_mag = sqrt(sq(lbm.u.x[n]) + sq(lbm.u.y[n]) + sq(lbm.u.z[n]));
				u_max_normal = fmax(u_max_normal, u_mag);
				
				if (x >= x_min && x <= x_max && y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
					if (u_mag > DANGER_VELOCITY) {
						danger_normal = true;
					}
				}
			}
		}

		print_info("Normal wind - Max velocity: " + to_string(u_max_normal, 4u) + " lu = " + to_string(u_max_normal * S_u, 1u) + " m/s");

#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
		//lbm.graphics.set_camera_free(float3(-0.464330f*(float)Nx, -0.791015f*(float)Ny,  2.138086f*(float)Nz), -117.1f,  24.4f,   8.0f);
		lbm.graphics.set_camera_free(float3(-0.464330f * (float)Nx, -0.791015f * (float)Ny, 2.138086f * (float)Nz), -117.1f, 24.4f, 8.0f);
		lbm.graphics.write_frame(get_exe_path() + "C/all_frames/"); // 始终保存
		if (danger_normal) {
			print_warning("DANGER ZONE DETECTED in view! Velocity exceeds 15 m/s!");
			lbm.graphics.write_frame(get_exe_path() + "破碎/danger_normal/"); // 危险时额外保存一份
		}
#endif

		print_info("=== CYCLE " + to_string(cycle) + ": Gust wind (1.5x) ===");

		// ========== 阵风 ==========
		lbm.u.read_from_device();
		for (ulong n = 0ull; n < lbm.get_N(); n++) {
			uint x = 0u, y = 0u, z = 0u;
			lbm.coordinates(n, x, y, z);
			if (x == 0u && lbm.flags[n] == TYPE_E) {
				lbm.u.x[n] = u_gust;
			}
		}
		lbm.u.write_to_device();
		lbm.run(500u, lbm_T);

		// ========== 速度统计（Gust wind）—— 仅监控视野内区域 ==========
		lbm.u.read_from_device();
		float u_max_gust = 0.0f;
		bool danger_gust = false;

		for (ulong n = 0ull; n < lbm.get_N(); n++) {
			if (lbm.flags[n] != TYPE_S) {
				uint x = 0u, y = 0u, z = 0u;
				lbm.coordinates(n, x, y, z);
				const float u_mag = sqrt(sq(lbm.u.x[n]) + sq(lbm.u.y[n]) + sq(lbm.u.z[n]));
				u_max_gust = fmax(u_max_gust, u_mag);
				
				if (x >= x_min && x <= x_max && y >= y_min && y <= y_max && z >= z_min && z <= z_max) {
					if (u_mag > DANGER_VELOCITY) {
						danger_gust = true;
					}
				}
			}
		}

		print_info("Gust wind - Max velocity: " + to_string(u_max_gust, 4u) + " lu = " + to_string(u_max_gust * S_u, 1u) + " m/s");

#if defined(GRAPHICS) && !defined(INTERACTIVE_GRAPHICS)
		//lbm.graphics.set_camera_free(float3(0.125100f * (float)Nx, -0.027440f * (float)Ny, -0.265927f * (float)Nz), -41.2f, 17.4f, 84.0f);
		lbm.graphics.set_camera_free(float3(-0.464330f * (float)Nx, -0.791015f * (float)Ny, 2.138086f * (float)Nz), -117.1f, 24.4f, 8.0f);
		lbm.graphics.write_frame(get_exe_path() + "C/all_frames/"); // 始终保存
		if (danger_gust) {
			print_warning("DANGER ZONE DETECTED in view! Velocity exceeds 15 m/s!");
			lbm.graphics.write_frame(get_exe_path() + "破碎/danger_gust/"); // 危险时额外保存一份
		}
#endif

		cycle++;
	}
} /**/