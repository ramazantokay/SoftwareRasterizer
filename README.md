# SoftwareRasterizer
Software Rasterizer - Forward Rendering Pipeline Implementation

This is a basic software rasterizer, implemented Modeling Transformation, Viewing Transformation, and Rasterization stages of the Forward Rendering Pipeline in C++.

Some of the properties of the Software Rasterizer as follows:
- Supported two types of projection transformations: orthographic projection and perspective projection.
- Used the midpoint algorithm to draw triangle edges and use the barycentric coordinates based algorithm to rasterize triangle faces.
- Liang-Barsky Clipping algorithm is used
- In both wireframe and solid modes, triangles whose backface is facing the viewer will be culled.
- Having strong OOP features
- Parsing XML scene files and rendering images in PPM format

### QuickStart

Download the repository

``` 
git clone https://github.com/ramazantokay/SoftwareRasterize.git
```
and run the following command on the terminal
```
make all
```
then, you can use the RayTracer with this command

```
./rasterize <sample_scenes/scene_file.xml>
```
It will render the scene and save it in its current directory with ppm format.

# Some provided scene 

**Horse and Mug**

![HorseMug](/assets/horse_and_mug_1.png)

**Filled Box**

![FilledBox](/assets/filled_box_4.png)

**Empty Box Culling Disabled**

![EmptyBox](/assets/empty_box_6.png)

**Empty Box Culling Enabled**

![EmptyBoxCE](/assets/empty_box_6_culling_enable.png)

**Empty Box Clipping**

![EmptyBoxClipping](/assets/empty_box_clipped_1_clipping.png)
