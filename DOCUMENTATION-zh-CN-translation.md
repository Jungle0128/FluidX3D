# FluidX3D 文档 - 如何入门？

## 0\. 安装 GPU 驱动程序和 OpenCL 运行时

（点击展开该部分）

*   **Windows**
    
    GPU
    
    *   下载并安装 AMD/Intel/Nvidia 显卡驱动程序，其中包含 OpenCL 运行时环境。
    *   重启。
    
    中央处理器
    
    *   下载并安装适用于 OpenCL 的英特尔 CPU 运行时（兼容 AMD/英特尔 CPU）。
    *   重启。
    
*   **Linux**
    
    AMD 显卡
    
    *   下载并安装 AMD GPU 驱动程序（内含 OpenCL 运行时），操作如下：
        
        ```bash
        sudo apt update && sudo apt upgrade -y
        sudo apt install -y g++ git make ocl-icd-libopencl1 ocl-icd-opencl-dev
        mkdir -p ~/amdgpu
        wget -P ~/amdgpu https://repo.radeon.com/amdgpu-install/6.4.2.1/ubuntu/noble/amdgpu-install_6.4.60402-1_all.deb
        sudo apt install -y ~/amdgpu/amdgpu-install*.deb
        sudo amdgpu-install -y --usecase=graphics,rocm,opencl --opencl=rocr
        sudo usermod -a -G render,video $(whoami)
        rm -r ~/amdgpu
        sudo shutdown -r now
        ```
    英特尔图形处理器
    
    *   自 Linux 内核 6.2 起，英特尔 GPU 驱动程序已预装，但其中不包含 OpenCL 运行时。
    *   OpenCL 运行时必须单独安装，使用以下命令：
        
        ```bash
        sudo apt update && sudo apt upgrade -y
        sudo apt install -y g++ git make ocl-icd-libopencl1 ocl-icd-opencl-dev intel-opencl-icd
        sudo usermod -a -G render $(whoami)
        sudo shutdown -r now
        ```
    英伟达图形处理器
    
    *   下载并安装 Nvidia GPU 驱动程序（内含 OpenCL 运行时），操作如下：
        
        ```bash
        sudo apt update && sudo apt upgrade -y
        sudo apt install -y g++ git make ocl-icd-libopencl1 ocl-icd-opencl-dev nvidia-driver-580
        sudo shutdown -r now
        ```
    中央处理器
    
    *   选项 1：使用以下命令下载并安装 oneAPI DPC++ 编译器和 oneTBB：
        
        ```bash
        export OCLV="oclcpuexp-2025.20.6.0.04_224945_rel"
        export TBBV="oneapi-tbb-2022.2.0"
        sudo apt update && sudo apt upgrade -y
        sudo apt install -y g++ git make ocl-icd-libopencl1 ocl-icd-opencl-dev
        sudo mkdir -p ~/cpurt /opt/intel/${OCLV} /etc/OpenCL/vendors /etc/ld.so.conf.d
        sudo wget -P ~/cpurt https://github.com/intel/llvm/releases/download/2025-WW27/${OCLV}.tar.gz
        sudo wget -P ~/cpurt https://github.com/uxlfoundation/oneTBB/releases/download/v2022.2.0/${TBBV}-lin.tgz
        sudo tar -zxvf ~/cpurt/${OCLV}.tar.gz -C /opt/intel/${OCLV}
        sudo tar -zxvf ~/cpurt/${TBBV}-lin.tgz -C /opt/intel
        echo /opt/intel/${OCLV}/x64/libintelocl.so | sudo tee /etc/OpenCL/vendors/intel_expcpu.icd
        echo /opt/intel/${OCLV}/x64 | sudo tee /etc/ld.so.conf.d/libintelopenclexp.conf
        sudo ln -sf /opt/intel/${TBBV}/lib/intel64/gcc4.8/libtbb.so /opt/intel/${OCLV}/x64
        sudo ln -sf /opt/intel/${TBBV}/lib/intel64/gcc4.8/libtbbmalloc.so /opt/intel/${OCLV}/x64
        sudo ln -sf /opt/intel/${TBBV}/lib/intel64/gcc4.8/libtbb.so.12 /opt/intel/${OCLV}/x64
        sudo ln -sf /opt/intel/${TBBV}/lib/intel64/gcc4.8/libtbbmalloc.so.2 /opt/intel/${OCLV}/x64
        sudo ldconfig -f /etc/ld.so.conf.d/libintelopenclexp.conf
        sudo rm -r ~/cpurt
        ```
        
    *   选项 2：使用以下命令下载并安装 PoCL：
        
        ```bash
        sudo apt update && sudo apt upgrade -y
        sudo apt install -y g++ git make ocl-icd-libopencl1 ocl-icd-opencl-dev pocl-opencl-icd
        ```
    
*   **安卓**
    
    ARM 显卡
    
    *   下载 Termux 并安装它。
    *   在 Termux 应用中运行：
        
        ```bash
        apt update && apt upgrade -y
        apt install -y clang git make
        ```

## 1\. 下载 FluidX3D

*   [下载](https://github.com/ProjectPhysX/FluidX3D/archive/refs/heads/master.zip)并解压源代码，或使用以下命令克隆：
    
    ```bash
    git clone https://github.com/ProjectPhysX/FluidX3D.git && cd FluidX3D
    ```
    
*   要更新 FluidX3D：
    *   备份您的更改。
    *   运行：
        
        ```bash
        git stash
        git pull origin master
        git stash pop
        ```
        

## 2\. 编译源代码

*   FluidX3D 没有"安装"过程，您需要自行编译源代码。
*   我已尽可能简化操作流程，本文档将引导您完成整个过程。不过，若具备 C++基础编程经验，将有助于您更好地理解配置脚本。
*   首先，按原样编译代码；这是标准的 FP32 基准测试用例。默认情况下，系统会自动选择最快的已安装 GPU。编译时间约为 5 秒。

### Windows

*   下载并安装 [Visual Studio Community](https://visualstudio.microsoft.com/de/vs/community/)。在 Visual Studio 安装程序中添加：
    *   桌面应用开发与 C++
    *   MSVC v142
    *   Windows 10 软件开发工具包
*   在 [Visual Studio Community](https://visualstudio.microsoft.com/de/vs/community/) 中打开 [`FluidX3D.sln`](FluidX3D.sln)。
*   点击 ► Local Windows Debugger 按钮进行编译和运行。
*   要选择特定的 GPU，请在 `FluidX3D` 文件夹中打开 Windows 命令提示符（在文件资源管理器的目录栏输入 `cmd` 并按 Enter ），然后运行 `bin\FluidX3D.exe 0` 以选择设备 `0`。若系统已[配置为多 GPU 模式](#the-lbm-class) ，也可通过 `bin\FluidX3D.exe 0 1 3 6` 同时选择多个 GPU。

### Linux / macOS / Android

*   使用以下命令编译并运行：
    
    ```bash
    chmod +x make.sh
    ./make.sh
    ```
    
*   编译需要支持 `C++17` 的 [`g++`](https://gcc.gnu.org/)（自 `8` 版本起支持，可通过 `g++ --version` 验证）。若已安装 [`make`](https://www.gnu.org/software/make/)（请通过 `make --version` 验证），编译将利用多核 CPU 加速；否则编译将降级为单核 CPU 运行。
*   要选择特定 GPU，请输入 `./make.sh 0` 进行编译+运行，或输入 `bin/FluidX3D 0` 在设备 `0` 上运行。若系统已[配置为多 GPU 模式](#the-lbm-class) ，也可通过 `bin/FluidX3D 0 1 3 6` 同时选择多个 GPU。
*   操作系统（Linux/macOS/Android）和 X11 支持（[`INTERACTIVE_GRAPHICS`](src/defines.hpp) 所需）将自动检测。若出现问题，您仍可手动在 [`make.sh`](make.sh#L13) 中选择 [`target=...`](make.sh#L13)。
*   在 macOS 和 Android 平台上，由于无法使用 X11，因此不支持 [`INTERACTIVE_GRAPHICS`](src/defines.hpp) 模式。不过您仍可使用 [`INTERACTIVE_GRAPHICS_ASCII`](src/defines.hpp) 模式，或通过常规 [`GRAPHICS`](src/defines.hpp) 模式将[视频渲染](#video-rendering)至硬盘。

## 3\. 浏览示例设置

*   现在打开 [`src/setup.cpp`](src/setup.cpp)。此处包含所有示例设置，每个示例都是用 C++编写的 `void main_setup() {...}` 函数块。取消其中一个的注释，可以从上到下依次尝试。
*   在 `main_setup()` 函数起始行中，显示"所需扩展项在 defines.hpp 中："，后接一列大写字母的扩展名。请前往 [`src/defines.hpp`](src/defines.hpp) 文件，将以下内容注释掉：
    
    ```c
    //#define BENCHMARK
    ```
    
    使用 `//` 注释掉所有扩展程序。然后，通过删除对应行前的 `//` 来取消注释所有设置所需的扩展程序。
*   最后，使用 ► Local Windows Debugger 按钮（Windows）或 `./make.sh`（Linux/macOS/Android） [编译](#2-compiling-the-source-code)并运行设置程序。
*   交互式图形窗口打开后，按键 P 可启动/暂停模拟，按键 H 可显示键盘控制与可视化设置的帮助菜单。
*   通过这种方式浏览一些示例配置，熟悉它们的代码结构并测试图形模式。

## 4\. [`INTERACTIVE_GRAPHICS`](src/defines.hpp)/[`_ASCII`](src/defines.hpp) 的键盘/鼠标控制

| 关键 | 功能 |
| --- | --- |
|  |  |
| P | 开始/暂停模拟 |
| H | 显示/隐藏键盘控制和可视化设置的帮助菜单 |
| EscAlt+F4 | 退出 |
|  |  |
| MouseIJ K L | 旋转摄像头 |
| Scrollwheel+ - | 缩放（中心摄像机模式）或摄像机移动速度（自由摄像机模式） |
| MouseclickU | 使用 Mouse 切换旋转模式，使用 I 、 J 、 K 、 L 实现角度吸附旋转 |
| F | 切换居中/自由摄像机模式 |
| WA S DSpace C | 移动自由摄像机 |
| Y X | 调整摄像机视场角 |
| R | 切换相机自动旋转 |
| G | 在控制台打印当前摄像机位置/旋转，作为可复制粘贴的命令 |
| V | 切换 VR 的立体渲染 |
| B | 切换 VR 眼镜/3D 电视模式以实现立体渲染 |
| N M | 调整立体渲染的视距 |
|  |  |
| 1 | 旗标线框模型/实体表面（若使用扩展功能，则包含实体单元上的力矢量或表面压力） |
| 2 | 速度场 |
| 3 | 流线型 |
| 4 | 涡度（速度着色的 Q 准则等值面） |
| 5 | 光栅化自由表面 |
| 6 | 光线追踪自由曲面 |
| 7 | 粒子 |
| T | 切换切片可视化模式 |
| Z | 切换字段可视化模式 |
| Q E | 在切片可视化模式下移动切片 |

## 5\. 编写自己的设置

### LBM 课程

*   初始化模拟框时，请使用调用
    
    ```c
    LBM lbm(Nx, Ny, Nz, nu, ...);
    ```
    
    构造函数。`Nx`/`Ny`/`Nz` 为网格分辨率，`nu` 为运动剪切粘度（单位为 [LBM 单位](#unit-conversion) ）。
*   要使用多个 GPU，请使用
    
    ```c
    LBM lbm(Nx, Ny, Nz, Dx, Dy, Dz, nu, ...);
    ```
    
    其中 `Dx`/`Dy`/`Dz` 分别表示每个空间方向上的域（GPU）数量。`Dx`×`Dy`×`Dz` 的乘积即为域（GPU）的总数。
*   只要 `lbm` 对象处于作用域内，您就可以访问内存。一旦它超出作用域，当前仿真相关的所有内存都会被重新释放。
*   网格分辨率 `Nx`/`Ny`/`Nz` 最终决定了 VRAM 的占用量。通常很难明确在何种分辨率下会突破 GPU 的 VRAM 容量上限。为此提供了一个辅助函数：
    
    ```c
    const uint3 lbm_N = resolution(float3(1.0f, 2.0f, 0.5f), 2000u);
    ```
    
    该函数以模拟盒所需的宽高比和以 MB 为单位的显存占用量为输入，返回网格分辨率作为 `uint3` 类型，包含 `.x`/`.y`/`.z` 三个分量。您也可以直接将 `uint3` 作为分辨率参数传入 LBM 构造函数：
    
    ```c
    LBM lbm(lbm_N, nu, ...);
    ```
    

### 单位换算

*   LBM 模拟采用与国际单位制不同的单位系统，其中密度 `ρ=1` 且速度 `u≈0.001-0.1`，因为浮点运算在接近 `1` 时精度最高。
*   为简化从国际单位制（SI）到长板英尺（LBM）单位的单位转换，提供了 [`units.hpp`](src/units.hpp) 结构体。通过调用
    
    ```c
    units.set_m_kg_s(lbm_length, lbm_velocity, lbm_density=1, si_length, si_velocity, si_density);
    ```
    
    基本单位转换因子\[m\]、\[kg\]、\[s\]经计算后存储于 `units` 结构体中。此后，可使用 [`src/units.hpp`](src/units.hpp) 中的任意转换函数在 SI 单位与 LBM 单位间进行双向转换，例如通过 `lbm_nu = units.nu(si_nu)` 将运动粘度从 SI 单位转换为 LBM 单位。
*   一个很好的入门示例就是" [奶牛的空气动力学](src/setup.cpp) "的设置。

### 初始条件与边界条件

*   若未显式设置，则默认所有单元格均采用默认值：`rho=1`、`u=0`、`flags=0`。
*   单个网格单元的初始/边界条件是在一个并行化循环中设置的，该循环遍历整个网格：
    
    ```c
    const uint Nx=lbm.get_Nx(), Ny=lbm.get_Ny(), Nz=lbm.get_Nz(); parallel_for(lbm.get_N(), [&](ulong n) { uint x=0u, y=0u, z=0u; lbm.coordinates(n, x, y, z);
    	// ...
    });
    ```
    
    在此循环中，可通过为 `lbm.rho[n]`、`lbm.u.x[n]`、`lbm.u.y[n]`、`lbm.u.z[n]` 及 `lbm.flags[n]` 赋值，分别设置每个网格单元的密度、速度及标志位。此处的 `n` 为线性化三维网格索引，通过函数 `lbm.coordinates(n, x, y, z)` 可对应转换为(`x`|`y`|`z`)坐标位置。
    *   例如，要在初始化循环中为模拟盒的四面设置实心边界，请编写：
        
        ```c
        if(x==0u||x==Nx-1u||y==0u||y==Ny-1u||z==0u||z==Nz-1u) {
        	lbm.flags[n] = TYPE_S;
        }
        ```
        
    *   对于基本几何体（如球体、椭球体、立方体、长方体、圆柱体、圆柱管、三角形、斜面或环面），可通过 [`shapes.hpp`](src/shapes.hpp) 中的函数进行设置。插入圆柱体的示例如下：
        
        ```c
        if(cylinder(x, y, z, lbm.center(), float3(axis_x, axis_z, axis_z), radius) {
        	lbm.flags[n] = TYPE_S;
        }
        ```
        
*   可用边界条件
    *   周期性边界
        *   所有未设置实心边界（`TYPE_S`）或其他边界类型的盒子侧面，均将保持周期性边界。
        *   若需严格质量守恒（例如线性管道中的流体流动），请采用周期性边界（即在模拟盒体侧面的网格单元上不设置任何边界类型），并通过体积力（相当于压力梯度）驱动流体流动。为此需启用（取消注释）[`VOLUME_FORCE`](src/defines.hpp) 扩展模块，并在 [LBM 构造器](#the-lbm-class)中设置体积力（`fx`|`fy`|`fz`）：
            
            ```c
            LBM lbm(Nx, Ny, Nz, nu, fx, fy, fz);
            ```
            
            这些单位体积内的力值不应超过 `0.001` 的数值。
    *   实体边界（`TYPE_S`）
        *   始终可用，无需进一步扩展。
        *   固体边界作为固体非移动物体，以非移动无滑移的中网格反弹边界形式实现。
            *   "静止" - 边界速度为0
            *   "无滑移反弹" - 边界处的流速为零（无滑移条件）
            *   "中网格" - 边界点精确位于边界单元与相邻流体单元的中间位置
        *   用法：
            
            ```c
            if(<condition defining the geometry>) {
            	lbm.flags[n] = TYPE_S;
            }
            ```
            
    *   平衡边界（指定密度和/或速度的 `TYPE_E`）
        *   您需要启用（取消注释）[`EQUILIBRIUM_BOUNDARIES`](src/defines.hpp) 扩展。
        *   平衡边界作为流入/流出边界发挥作用。`lbm.flags[n] = TYPE_E` 单元强制执行初始化时指定的密度 `lbm.rho[n]` 和/或速度 `lbm.u.x[n]`/`lbm.u.y[n]`/`lbm.u.z[n]`，并吸收任何传入的激波。若未指定密度/速度值，则强制采用默认值（`rho=1`，`u=0`）。
        *   用法：对于特定的流入/流出单元，设置标志 `lbm.flags[n] = TYPE_E`，并在相同单元上指定密度 `lbm.rho[n]` 不等于 `1` 或速度 `lbm.u.x[n]`/`lbm.u.y[n]`/`lbm.u.z[n]` 不同于 `0`，或两者组合。
            
            ```c
            if(<condition defining the geometry>) {
            	lbm.flags[n] = TYPE_E;
            	lbm.rho[n] = <...>; // TYPE_E cells will retain this fixed density
            	lbm.u.x[n] = <...>; // TYPE_E cells will retain this fixed velocity
            	lbm.u.y[n] = <...>; // TYPE_E cells will retain this fixed velocity
            	lbm.u.z[n] = <...>; // TYPE_E cells will retain this fixed velocity
            }
            ```
            
    *   移动固体边界（`TYPE_S` 且速度不为零）
        *   您需要启用（取消注释）[`MOVING_BOUNDARIES`](src/defines.hpp) 扩展。
        *   移动固体边界既可作为固体移动物体，也可作为具有固定流速的流入/流出装置（类似于液压活塞）。它们会反射任何入射冲击波。
        *   移动固体边界被实现为移动无滑移的中网格反弹边界。
            *   "移动" - 边界速度为指定的非零值（速度的狄利克雷边界条件）。
            *   "无滑移反弹" - 边界处的流速直接等于非零边界速度（无滑移条件）
            *   "中网格" - 边界点精确位于边界单元与相邻流体单元的中间位置
        *   用法：对于特定的移动固体单元，设置标志 `lbm.flags[n] = TYPE_S`，并在同一单元上指定非零速度 `lbm.u.x[n]`/`lbm.u.y[n]`/`lbm.u.z[n]`。
            
            ```c
            if(<condition defining the geometry>) {
            	lbm.flags[n] = TYPE_S;
            	lbm.u.x[n] = <...>; // TYPE_S cells will retain this fixed velocity
            	lbm.u.y[n] = <...>; // TYPE_S cells will retain this fixed velocity
            	lbm.u.z[n] = <...>; // TYPE_S cells will retain this fixed velocity
            }
            ```
            
        *   您还可通过 [`SURFACE`](#surface-extension) 扩展模块，将移动固体边界作为固定流速的流入源。此时除在 `TYPE_S` 单元格中指定流速外，还需使用 `TYPE_F` 标记初始化毗邻 `TYPE_S` 单元格的单元层。

### 运行模拟

*   调用 `lbm.run()`（不带输入参数时表示无限时间步长）来初始化并执行设置，或调用 `lbm.run(time_steps)` 仅执行特定数量的时间步长。
*   若您拥有一个[更复杂的模拟循环](#video-rendering) ，其中需要周期性计算时间步长并渲染视频图像或导出数据，请务必在该循环前添加 `lbm.run(0u)`。此操作将初始/边界条件从 CPU 内存复制至 GPU 显存，并在 GPU 上初始化模拟，同时跳过时间步长的计算。若未进行初始化，VRAM 中将尚无可用于渲染的数据。

### 加载 .stl 文件

*   对于更复杂的几何体，您可以加载 `.stl` 三角网格，并在 GPU 上将其体素化到笛卡尔模拟网格中。
*   在 `FluidX3D/src/` 文件夹旁创建一个 `FluidX3D/stl/` 文件夹，并从 [Thingiverse](https://www.thingiverse.com/) 等网站下载几何模型，或自行创建。
*   仅支持二进制 `.stl` 文件格式。若需转换其他格式或拆分复合几何体（如直升机机身与旋翼），建议在 Windows 系统使用 [Microsoft 3D Builder](https://apps.microsoft.com/store/detail/3d-builder/9WZDNCRFJ3T6)，或在 Windows/Linux 系统使用 [Blender](https://www.blender.org/)。
*   直接加载并体素化简单的 `.stl` 文件
    
    ```c
    lbm.voxelize_stl(get_exe_path()+"../stl/mesh.stl", center, rotation, size);
    ```
    
    这将自动将网格重新定位/缩放至指定中心。使用 `lbm.center()` 对准模拟盒中心，或通过 `+float3(offset_x, offset_y, offset_z)` 添加偏移量。可按如下方式生成并乘积旋转矩阵（示例：先绕 z 轴旋转 180°，再绕 x 轴旋转 90°）：
    
    ```c
    float3x3 rotation = float3x3(float3(1, 0, 0), radians(90.0f))*float3x3(float3(0, 0, 1), radians(180.0f));
    ```
    
*   要加载包含多个部件的复合几何体，且不进行自动网格重新定位/缩放，请使用
    
    ```c
    Mesh* mesh_1 = read_stl(const string& path, const float scale=1.0f, const float3x3& rotation=float3x3(1.0f), const float3& offset=float3(0.0f)); // load mesh without automatic repositioning/rescaling
    Mesh* mesh_2 = read_stl(const string& path, const float scale=1.0f, const float3x3& rotation=float3x3(1.0f), const float3& offset=float3(0.0f));
    mesh_1->scale(const float scale); // manually scale meshes
    mesh_2->scale(const float scale);
    mesh_1->translate(const float3& translation); // manually reposition meshes
    mesh_2->translate(const float3& translation);
    lbm.voxelize_mesh_on_device(mesh_1); // voxelize meshes on GPU
    lbm.voxelize_mesh_on_device(mesh_2);
    
    ```
    
    从 `.stl` 文件加载网格，同时手动缩放/重新定位网格的所有部分，最后在 GPU 上将其体素化。
*   为辅助网格重新定位，提供以下函数：- \`lbm.center()\` 获取模拟盒中心坐标- \`mesh->pmin\`/\`mesh->pmax\` 获取网格最小/最大边界框坐标- \`mesh->pmax\` 获取网格最大边界框坐标每个坐标均为 \`float3\` 类型，包含 \`(x|y|z)\` 三个分量。3Epmax，均为 `float3` 类型，包含(x|y|z)三个分量。
*   旋转几何体需要定期重新体素化，大约每 1-10 个 LBM 时间步长进行一次。在 [`main_setup()`](src/setup.cpp) 函数的主模拟循环中，首先旋转三角网格，然后在 GPU 上重新体素化，最后计算若干 LBM 时间步长：
    
    ```c
    const uint lbm_T = 100000u; // number of LBM time steps to simulate
    const uint lbm_dt = 4u; // number of LBM time steps between each mesh revoxelization
    lbm.run(0u); // initialize simulation
    while(lbm.get_t()<lbm_T) { // main simulation loop
    	mesh->rotate(float3x3(float3(0, 0, 1), lbm_omega*(float)lbm_dt)); // rotate the triangle mesh
    	lbm.voxelize_mesh_on_device(mesh, TYPE_S, center, float3(0.0f), float3(0.0f, 0.0f, lbm_omega)); // revoxelize the rotated triangle mesh, provide the instantaneous angular velocity vector for moving boundaries
    	lbm.run(lbm_dt); // run lbm_dt LBM time steps
    }
    ```
    
    此处 `lbm_omega` 表示每时间步的角速度（单位：弧度），`lbm_dt` 表示两次重采样之间的模拟时间步数，而 `float3(0.0f, 0.0f, lbm_omega)` 表示沿旋转轴方向的瞬时角速度向量。外层单元格在两次旋转重采样间隔内的最大位移不应超过 `1` 个单元格；据此将 `lbm_omega = lbm_u/lbm_radius` 进行相应设置。
*   请查看" [塞斯纳 172](src/setup.cpp)"和" [贝尔 222](src/setup.cpp)"的配置示例。

### 视频渲染

*   对于视频渲染，请禁用（注释掉）[`INTERACTIVE_GRAPHICS`](src/defines.hpp) 和 [`INTERACTIVE_GRAPHICS_ASCII`](src/defines.hpp)，并启用 [`GRAPHICS`](src/defines.hpp)。
*   将视频分辨率设置为 [`GRAPHICS_FRAME_WIDTH`](src/defines.hpp)/[`GRAPHICS_FRAME_HEIGHT`](src/defines.hpp)，并将背景色设置为 [`GRAPHICS_BACKGROUND_COLOR`](src/defines.hpp)。您还可以在那里调整其他 [`GRAPHICS_...`](src/defines.hpp) 选项，例如半透明渲染模式，或通过 [`GRAPHICS_U_MAX`](src/defines.hpp) 调整速度的颜色标尺。
*   在 [`main_setup()`](src/setup.cpp) 函数中，用于渲染视频的基本循环如下所示：
    
    ```c
    lbm.graphics.visualization_modes = VIS_FLAG_LATTICE|VIS_Q_CRITERION; // set visualization modes, see all available visualization mode macros (VIZ_...) in defines.hpp
    const uint lbm_T = 10000u; // number of LBM time steps to simulate
    lbm.run(0u, lbm_T); // initialize simulation
    while(lbm.get_t()<lbm_T) { // main simulation loop
    	if(lbm.graphics.next_frame(lbm_T, 25.0f)) { // render enough frames for 25 seconds of 60fps video
    		lbm.graphics.set_camera_free(float3(2.5f*(float)Nx, 0.0f*(float)Ny, 0.0f*(float)Nz), 0.0f, 0.0f, 50.0f); // set camera to position 1
    		lbm.graphics.write_frame(get_exe_path()+"export/camera_1/"); // export image from camera position 1
    		lbm.graphics.set_camera_centered(-40.0f, 20.0f, 78.0f, 1.25f); // set camera to position 2
    		lbm.graphics.write_frame(get_exe_path()+"export/camera_2/"); // export image from camera position 2
    	}
    	lbm.run(1u, lbm_T); // run 1 LBM time step
    }
    ```
    
*   为确定合适的摄像机位置，请以低分辨率在 [`INTERACTIVE_GRAPHICS`](src/defines.hpp) 模式下运行模拟，旋转/移动摄像机至目标位置，点击 Mouse 禁用鼠标旋转，按下 G 将当前摄像机设置以可复制粘贴的命令形式输出至控制台。将 Alt + Tab 输出至控制台，通过鼠标选中并右键点击复制摄像机定位命令，随后将其粘贴至 [`main_setup()`](src/setup.cpp) 函数中。
*   要使摄像机沿着一组预设的关键帧位置平滑移动，请使用 `catmull_rom` 样条曲线：
    
    ```c
    while(lbm.get_t()<=lbm_T) { // main simulation loop
    	if(lbm.graphics.next_frame(lbm_T, 30.0f)) {
    		const float t = (float)lbm.get_t()/(float)lbm_T;
    		vector<float3> camera_positions = {
    			float3(-0.282220f*(float)Nx,  0.529221f*(float)Ny,  0.304399f*(float)Nz),
    			float3( 0.806921f*(float)Nx,  0.239912f*(float)Ny,  0.436880f*(float)Nz),
    			float3( 1.129724f*(float)Nx, -0.130721f*(float)Ny,  0.352759f*(float)Nz),
    			float3( 0.595601f*(float)Nx, -0.504690f*(float)Ny,  0.203096f*(float)Nz),
    			float3(-0.056776f*(float)Nx, -0.591919f*(float)Ny, -0.416467f*(float)Nz)
    		};
    		vector<float> camera_rx = {
    			 116.0f,
    			  25.4f,
    			 -10.6f,
    			 -45.6f,
    			 -94.6f
    		};
    		vector<float> camera_ry = {
    			  26.0f,
    			  33.3f,
    			  20.3f,
    			  25.3f,
    			 -16.7f
    		};
    		const float camera_fov = 90.0f;
    		lbm.graphics.set_camera_free(catmull_rom(camera_positions, t), catmull_rom(camera_rx, t), catmull_rom(camera_ry, t), camera_fov);
    		lbm.graphics.write_frame(get_exe_path()+"export/");
    	}
    	lbm.run(1u, lbm_T);
    }
    ```
    
*   可通过 [`VIS_...`](src/defines.hpp) 宏指定可视化模式，格式为 `lbm.graphics.visualization_modes` 。您还可设置 `lbm.graphics.slice_mode`（`0`\=不切片，`1`\=x 轴，`2`\=y 轴，`3`\=z 轴，`4`\=xz 轴，`5`\=xyz 轴，`6`\=yz 轴，`7`\=xy）并通过 `lbm.graphics.slice_x`/`lbm.graphics.slice_y`/`lbm.graphics.slice_z` 重新定位切片。
*   导出的帧将在其名称中自动分配当前模拟时间步长，格式为 `bin/export/image-123456789.png` 。
*   要将渲染的 `.png` 图像转换为视频，请使用 [FFmpeg](https://ffmpeg.org/)：
    
    ```bash
    ffmpeg -framerate 60 -pattern_type glob -i "export/*/image-*.png" -c:v libx264 -pix_fmt yuv420p -b:v 24M "video.mp4"
    ```
    

### 数据导出

*   在任何时刻，您都可以通过以下方式将体积数据导出为二进制 `.vtk` 文件：
    
    ```c
    lbm.rho.write_device_to_vtk(); // density
    lbm.u.write_device_to_vtk(); // velocity
    lbm.flags.write_device_to_vtk(); // flags
    lbm.F.write_device_to_vtk(); // force, only for FORCE_FIELD extension
    lbm.phi.write_device_to_vtk(); // fill fraction, only for SURFACE extension
    lbm.T.write_device_to_vtk(); // temperature, only for TEMPERATURE extension
    lbm.write_mesh_to_vtk(const Mesh* mesh); // for exporting triangle meshes
    ```
    
*   这些函数首先将数据从 GPU(s)中提取到 CPU 内存中，然后将其写入硬盘。
*   若通过 `units.set_m_kg_s(...)` 指定了[单位转换](#unit-conversion) ，导出的 `.vtk` 文件中的数据将自动转换为国际单位制（SI）单位。
*   导出的文件将自动在其名称中添加当前模拟时间步长，格式为 `bin/export/u-123456789.vtk`。
*   请注意，这些体积文件的文件大小可能非常庞大，单个文件可能达到数十 GB。
*   您可以在 [ParaView](https://www.paraview.org/) 等软件中查看/评估 `.vtk` 文件。
*   建议直接在 [`main_setup()`](src/setup.cpp) 函数中使用 C++功能提取所需数据，并选择性地仅将该数据写入硬盘。因此，调用 `lbm.u.read_from_device()` 将数据从 GPU 复制到 CPU 内存，之后即可直接访问，例如：
    
    ```c
    const float lbm_velocity_x = lbm.u.x[lbm.index(x, y, z)];
    ```
    
    获取位置（`x`|`y`|`z`）处的 x 方向速度，单位为 [LBM 单位](#unit-conversion) 。
*   要将 [LBM 单位的流速转换为 SI 单位](#unit-conversion) ，请使用
    
    ```c
    const float si_velocity_x = units.si_u(lbm_velocity_x);
    ```
    
    在使用 `units.set_m_kg_s(...)` 完成 [单位转换](#unit-conversion) 之后。
*   您还可以通过以下方式将 `.stl` 三角网格导出为二进制 `.vtk` 文件：
    
    ```c
    lbm.write_mesh_to_vtk(const Mesh* mesh);
    ```
    

### 升力/阻力

*   启用（取消注释）[`FORCE_FIELD`](src/defines.hpp) 扩展。该扩展允许对每个实体单元（`TYPE_S`）单独计算边界力，同时可在每个流体单元上施加独立体积力（此处未使用）。
*   在 [`main_setup()`](src/setup.cpp) 函数中，使用独特的标志组合对网格进行体素化处理，例如 `(TYPE_S|TYPE_X)` 或 `(TYPE_S|TYPE_Y)` 或 `(TYPE_S|TYPE_X|TYPE_Y)`，以此区别于可能用于定义其他几何体的所有 `(TYPE_S)` 单元，并计算其质心：
    
    ```c
    lbm.voxelize_mesh_on_device(mesh, TYPE_S|TYPE_X); // voxelize mesh with unique flag combination
    const float3 lbm_com = lbm.object_center_of_mass(TYPE_S|TYPE_X); // object center of mass in LBM unit coordinates
    ```
    
*   在 [`main_setup()`](src/setup.cpp) 函数的主模拟循环调用中，对属于该对象的所有边界单元进行求和：
    
    ```c
    const float3 lbm_force = lbm.object_force(TYPE_S|TYPE_X); // force on object
    const float3 lbm_torque = lbm.object_torque(lbm_com, TYPE_S|TYPE_X); // torque on object around lbm_com rotation point
    ```
    
    这些函数对属于该对象的所有标记为 `(TYPE_S|TYPE_X)` 的单元进行求和。求和操作在显存中通过 GPU 加速完成，最终仅将结果复制到 CPU 内存。
*   您还可访问单个网格单元的力场。请注意，将整个 `lbm.F` 力场从 GPU 显存复制到 CPU 内存会很耗时：
    
    ```c
    lbm.F.read_from_device(); // copy entire force field from GPU VRAM to CPU RAM (slow)
    lbm_force_x_n = lbm.F.x[lbm.index(x, y, z)]; // access force at one particular grid cell with integer coordinates x, y, z
    ```
    
*   最后， [将英制单位转换为国际单位制](#unit-conversion) ，使用
    
    ```c
    const float si_force_x = units.si_F(lbm_force.x);
    ```
    
    在使用 `units.set_m_kg_s(...)` 完成 [单位转换](#unit-conversion) 之后。
*   请参阅"[Ahmed 模型](src/setup.cpp) "的示例设置。需注意在高度湍流工况下，计算得到的体积力可能过大（最大可达 2 倍），因为即使采用高分辨率仍无法完全捕捉湍流边界层。此时仍需采用壁面函数。

## 6\. LBM 的进一步扩展

至此，您已通过扩展功能熟悉了[额外的边界类型](#initial-and-boundary-conditions) ，包括[`体积力`](src/defines.hpp) 、[`FORCE_FIELD`](src/defines.hpp)、[`EQUILIBRIUM_BOUNDARIES`](src/defines.hpp) 以及 [`MOVING_BOUNDARIES`](src/defines.hpp)。其余可用模型扩展功能简要概述如下：

### [`SURFACE`](src/defines.hpp) 扩展

*   要模拟自由水面，请启用（取消注释）[`SURFACE`](src/defines.hpp) 扩展。
*   所有单元格随后获得 3 个额外标记：`TYPE_F`（流体）、`TYPE_I`（界面）和 `TYPE_G`（气体）。流体单元采用常规 LBM 方法计算。界面单元需考虑额外表面张力作用——当 [LBM 构造函数](#the-lbm-class)中表面张力系数 `sigma` 设置为大于 `0` 时；界面厚度始终为 1 个单元层。气体单元完全不参与模拟，实质上被视为真空环境。
*   若未在[初始条件](#initial-and-boundary-conditions)中另行设置，所有单元格默认初始化为 `TYPE_G`。作为初始条件，将所有应为流体状态的单元格设置为
    
    ```c
    lbm.flags[n] = TYPE_F;
    ```
    
    接口层将在初始化过程中通过 `lbm.run(0u)` 自动初始化。
*   除了 3 个标志位外，每个网格单元还会被分配一个填充水平 `lbm.phi[n]`：流体单元（`TYPE_F`）填充为 `1`，界面单元（`TYPE_I`）填充为 `]0,1[`，气体单元（`TYPE_G`）填充为 `0`。除单元标志外，您还可在初始化时设置此填充级别。切勿忘记设置单元标志。若未手动设置 `lbm.phi[n]`，系统将自动初始化：所有流体单元设为 `phi=1`，所有界面单元设为 `phi=0.5`，所有气体单元设为 `phi=0`。
*   对于流体流入边界条件，请使用[移动固体边界](#initial-and-boundary-conditions) ，并在与 `TYPE_S` 单元相邻且具有非零速度的单元中，额外设置 `TYPE_F` 标记。
*   参见" [水力跃流](src/setup.cpp) "或" [雨滴冲击](src/setup.cpp) "示例装置。

### [`温度`](src/defines.hpp) 扩展

*   通过 [`TEMPERATURE`](src/defines.hpp) 扩展，FluidX3D 能够模拟热对流流动。该扩展同时自动启用 [`VOLUME_FORCE`](src/defines.hpp) 扩展。
*   在 [LBM 构造函数](#the-lbm-class)中，您需要设置体积力（`fx`|`fy`|`fz`）、热扩散系数 `alpha` 以及热膨胀系数 `beta`，所有参数均需[采用 LBM 单位](#unit-conversion) ：
    
    ```c
    LBM lbm(Nx, Ny, Nz, nu, fx, fy, fz, 0.0f, alpha, beta); // the "0.0f" is for the surface tension coefficient sigma which is not used here and has to remain 0
    ```
    
*   通过扩展，每个网格单元将额外分配一个温度值 `lbm.T[n]`（ [以 LBM 单位表示](#unit-conversion) ）。默认温度值为 `1`（LBM 单位）。
*   要设置温度边界条件，请使用标志 `TYPE_T`，并为相同单元格赋予不等于 `1` 的温度：
    
    ```c
    lbm.flags[n] = TYPE_T; // make the cell n a temperature boundary
    lbm.T[n] = 1.2f; // set this temperature boundary hotter than average
    ```
    
*   请参见" [雷利-贝纳德对流](src/setup.cpp) "和" [热对流](src/setup.cpp) "两个实例设置。

### [`SUBGRID`](src/defines.hpp) 扩展

*   流体流动的特征由雷诺数决定。
    
    *Re* = *x*·*u*∕*nu*
    
    具有特征长度尺度 `x`、特征速度 `u` 和运动剪切粘度 `nu`。较大的长度尺度、较大的速度或较小的粘度均意味着更大的雷诺数。
*   雷诺数是一个无量纲数。当数值较低（*Re* < 2300）时，表示层流；数值较高（*Re* > 2900）时，表示湍流。其间为过渡流场。
*   当雷诺数 *Re* 远大于 100000 时，LBM 求解器会变得[不稳定](#7-suitable-parameters-and-simulation-instability) 。这是因为流场中可能存在微小且高速旋转的涡旋，而过高的流速和剪切率会导致模拟结果失效。
*   为解决此问题，可采用亚网格模型来模拟小于单个网格单元的涡流。其原理是在剪切率较大且假定存在大量微小涡流的区域增加有效粘度。巧合的是，高剪切率与低粘度区域往往引发不稳定现象，因此通过提升该区域的有效粘度可维持模拟的稳定性。
*   FLuidX3D 中的次网格模型采用 Smagorinsky-Lilly 模型。可通过 [`SUBGRID`](src/defines.hpp) 扩展启用该功能。
*   此扩展不会带来额外的性能开销。

### [`粒子`](src/defines.hpp) 扩展

*   默认情况下，LBM 是一种基于网格的模拟，因此不存在粒子。
*   但 [`PARTICLES`](src/defines.hpp) 扩展允许向模拟中添加粒子，既可作为被动示踪剂，也可作为可进行浮选/沉降的双向耦合粒子。
*   对于被动追踪器，仅需启用 [`PARTICLES`](src/defines.hpp) 扩展，并在 [LBM 构造函数](#the-lbm-class)中直接添加粒子数量：
    
    ```c
    LBM lbm(Nx, Ny, Nz, nu, 50000u); // this will create 50000 particles
    ```
    
*   然后，在[初始化](#initial-and-boundary-conditions)阶段，遍历所有粒子（在遍历所有网格单元的初始化循环之外）：
    
    ```c
    uint seed = 42u;
    for(ulong n=0ull; n<lbm.particles->length(); n++) {
    	lbm.particles->x[n] = random_symmetric(seed, 0.5f*lbm.size().x); // this will palce the particles randomly anywhere in the simulation box
    	lbm.particles->y[n] = random_symmetric(seed, 0.5f*lbm.size().y);
    	lbm.particles->z[n] = random_symmetric(seed, 0.5f*lbm.size().z);
    }
    ```
    
*   请注意，粒子所在的位置（`0`|`0`|`0`）对应于模拟盒的中心。
*   对于双向耦合粒子，需额外启用 [`VOLUME_FORCE`](src/defines.hpp) 和 [`FORCE_FIELD`](src/defines.hpp) 扩展，并在 [LBM 构造函数](#the-lbm-class) 中添加粒子密度（ [以 LBM 单位表示](#unit-conversion) ）且不等于 `1`：
    
    ```c
    LBM lbm(Nx, Ny, Nz, nu, 50000u, 1.2f); // this will create 50000 particles that are more dense than the fluid and will sink to the bottom
    ```
    

## 7\. 适宜参数与模拟不稳定性

*   在速度场或流线可视化过程中，有时会出现模糊现象，或是类似快速生长的白色晶体从某一点膨胀开来，最终填满整个模拟区域。这便是数值不稳定现象——即当速度值变为 `NaN` 或 `Inf` 时所致。
*   通常，不稳定性的根源在于参数选择不当：
    *   密度过高/过低（理想情况下应始终非常接近1）
    *   速度过高 `u`（在整个区域内绝不能超过 `0.57`，理想值应在 `0.075` 左右，但可低至 `0.001`）
    *   运动剪切粘度 `nu` 过低（理想值接近 `1/6`，当其非常接近 `0` 时会变得不稳定（此时需启用 [`SUBGRID`](src/defines.hpp) 扩展），且不应超过 `3`）
    *   单位体积应力过高（`fx`|`fy`|`fz`）（其大小不应超过 `0.001`）
    *   表面张力系数 `sigma` 过高（不应超过 `0.1`）
*   对于 LBM 模拟而言，最佳参数化设置本身就是一门艺术，需要经过实践才能掌握。