// MARCHING CUBES RENDERING LIBRARY
// WRITTEN BY ANDRE SIER, 2010, ADAPTED BY SIMON EVER-HALE TO NOT REQUIRE EXTERNAL LIBRARIES AND WORK WITH PROCESSING.JS
// IT DOESN'T WORK GREAT THO :(
public class GRIDCELL {
  public XYZ p[]=new XYZ[8];
  public float val[] = new float[8];
  public GRIDCELL() {
    for (int i=0; i<8;i++) {
      p[i] = new XYZ(); 
      val[i] = 0.0f;
    }
  }
}



/*
 marching cubes implementation by andre sier in processing, 2010
   http://s373.net/code/marchingcubes
 adapted code from paul bourke's polygonizing a scale field (marching cubes), 
   http://paulbourke.net/geometry/polygonise/
 stl export code adapted from marius watz, 
   http://workshop.evolutionzone.com/unlekkerlib/
 */
//import java.io.*;
//import java.net.URISyntaxException;
//import java.nio.ByteBuffer;
//import java.nio.ByteOrder;
//import java.nio.FloatBuffer;

import processing.core.PApplet;
import processing.core.PConstants;

public class MarchingCubes {
  public PApplet applet;
  String VERSION = "s373.marchingcubes - 0.1.0 ";
  final String VERSIONURL = " - http://s373.net/code/marchingcubes \n";

  // obj data
  public GRIDCELL grid;
  public float isolevel;
  public TRIANGLE[] triangles;

  public int gx, gy, gz, numxyz, gxgy;
  public float data[];
  public float themin, themax;
  public ArrayList<TRIANGLE> trilist;
  public int ntri = 0;
  public boolean invertnormals = true;
  public boolean closesides = true;

  // draw data
  public XYZ worlddim, worldstride, worldcenter, datastride;
  public String info;

  //// file data
  //public String filename;
  //public File file;
  //public byte[] header, byte4;
  //public ByteBuffer buf;


  public MarchingCubes(PApplet pa, float sx, float sy, float sz, int x,
      int y, int z) {
    applet = pa;
    initResolution(x, y, z);
    setWorldDim(sx, sy, sz);

    triangles = new TRIANGLE[5];
    for (int i = 0; i < triangles.length; i++) {
      triangles[i] = new TRIANGLE();
    }
    isolevel = 0.0025f;
    grid = new GRIDCELL();
    polygoniseData();

  }

  public void setWorldDim(float x, float y, float z) {
    worlddim = new XYZ(x, y, z);
    worldstride = new XYZ(gx / x, gy / y, gz / z);
    datastride = new XYZ(x / gx, y / gy, z / gz);
    worldcenter = new XYZ(x/2.0f, y/2.0f, z/2.0f);
    //System.out.print("\nworlddim: " + x + " " + y + " " + z + "\n");
    //System.out.print("worldcenter: " + worldcenter.x + " " + worldcenter.y
    //    + " " + worldcenter.z + "\n");
    //System.out.print("worldstride: " + worldstride.x + " " + worldstride.y
    //    + " " + worldstride.z + "\n");
    //System.out.print("datastride: " + datastride.x + " " + datastride.y
    //    + " " + datastride.z + "\n");
  }

  public void initResolution(int x, int y, int z) {
    gx = x;
    gy = y;
    gz = z;
    gxgy = x * y;
    numxyz = x * y * z;
    data = new float[numxyz];
    for (int i = 0; i < numxyz; i++) {
      data[i] = random(0.0f, 0.7f);
    }
    trilist = new ArrayList<TRIANGLE>();
  }

  public String getinfo() {
    info = "tris: " + ntri + " volume: " + gx + " " + gy
    + " " + gz + " cells: " + numxyz + " iso: " + isolevel;
    return info;
  }

  public XYZ getCenter(){
    return worldcenter;
  }
  
  
  /*
   * ADDING VALUES methods
   */
  
  public final int getIndex(XYZ ptpos){
    int cx = (int) (ptpos.x * worldstride.x);
    int cy = (int) (ptpos.y * worldstride.y);
    int cz = (int) (ptpos.z * worldstride.z);
    if (closesides) {
      if (cx < 1) {
        cx = 1;
      }
      if (cy < 1) {
        cy = 1;
      }
      if (cz < 1) {
        cz = 1;
      }
      if (cx >= gx - 1) {
        cx = gx - 2;
      }
      if (cy >= gy - 1) {
        cy = gy - 2;
      }
      if (cz >= gz - 1) {
        cz = gz - 2;
      }
    } else {
      if (cx < 0) {
        cx = 0;
      }
      if (cy < 0) {
        cy = 0;
      }
      if (cz < 0) {
        cz = 0;
      }
      if (cx >= gx) {
        cx = gx - 1;
      }
      if (cy >= gy) {
        cy = gy - 1;
      }
      if (cz >= gz) {
        cz = gz - 1;
      }
    }
    return cx + cy * gx + cz * gx * gy;
  }
  
  
  public void addIsoPoint(float val, XYZ ptpos) {
    final int idx = getIndex(ptpos);
    data[idx] += val;
  }

  public void addCube(float val, XYZ ptpos, XYZ ptdim) {
    // replace this getindex
    
    int cx = (int) (ptpos.x * worldstride.x); // // width * gx);
    int cy = (int) (ptpos.y * worldstride.y);
    int cz = (int) (ptpos.z * worldstride.z);

    if (closesides) {
      if (cx < 1) {
        cx = 1;
      }
      if (cy < 1) {
        cy = 1;
      }
      if (cz < 1) {
        cz = 1;
      }
      if (cx >= gx - 1) {
        cx = gx - 2;
      }
      if (cy >= gy - 1) {
        cy = gy - 2;
      }
      if (cz >= gz - 1) {
        cz = gz - 2;
      }
    } else {
      if (cx < 0) {
        cx = 0;
      }
      if (cy < 0) {
        cy = 0;
      }
      if (cz < 0) {
        cz = 0;
      }
      if (cx >= gx) {
        cx = gx - 1;
      }
      if (cy >= gy) {
        cy = gy - 1;
      }
      if (cz >= gz) {
        cz = gz - 1;
      }
    }

    int dimx = (int) (ptdim.x * worldstride.x); // // width * gx);
    int dimy = (int) (ptdim.y * worldstride.y);
    int dimz = (int) (ptdim.z * worldstride.z);

    if (dimx < 2) {
      dimx = 2;
    }
    if (dimy < 2) {
      dimy = 2;
    }
    if (dimz < 2) {
      dimz = 2;
    }
    addCube(val, cx, cy, cz, dimx, dimy, dimz);
  }

  private void addCube(float val, int centerx, int centery, int centerz,
      int dimx, int dimy, int dimz) {

    final int hx = dimx / 2;
    final int hy = dimy / 2;
    final int hz = dimz / 2;

    int sx = centerx - hx;
    int sy = centery - hy;
    int sz = centerz - hz;

    if (sx < 1) {
      sx = 1;
    }
    if (sy < 1) {
      sy = 1;
    }
    if (sz < 1) {
      sz = 1;
    }

    final int tx = Math.min(centerx + hx, gx - 1);
    final int ty = Math.min(centery + hy, gy - 1);
    final int tz = Math.min(centerz + hz, gz - 1);

    // println("adding cube: "+val+" "+sz+" "+tz+" "+dimz);

    int idx = 0;
    for (int i = sx; i < tx; i++) {
      for (int j = sy; j < ty; j++) {
        for (int k = sz; k < tz; k++) {
          idx = i + j * gx + k * gx * gy;
          data[idx] += val;
        }
      }
    }
  }

  public void zeroData() {
    themin = 0;
    themax = 0;
    for (int i = 0; i < numxyz; i++) {
      data[i] = 0.0f;
    }
  }

  public void setRndData(float mn, float mx) {
    for (int i = 0; i < numxyz; i++) {
      data[i] = random(mn, mx);
    }
  }
  public void setRndData(float mx) {
    for (int i = 0; i < numxyz; i++) {
      data[i] = random(mx);
    }
  }

  public void setData(float d[]) {
    for (int i = 0; i < numxyz; i++) {
      data[i] = d[i];
    }
  }

  public void addData(float d[]) {
////phps add close sides here
    for (int i = 0; i < numxyz; i++) {
      data[i] += d[i];
    }
  }

  public float[] getData() {
    return data;
  }

  public boolean isEmpty(){
    boolean empty = true;
    for(int i=0; i< numxyz; i++){
      if(data[i]>0){
        empty = false;
        break;
      }
    }
    return empty;
  }

  public void multData(final float v){
    for (int i = 0; i < numxyz; i++) {
      data[i] *= v;
    }
  }
  
  public void normalizeDataTo(final float v){
    checkMinMax();
    float dist = themax - themin;
    final float dst = v / dist;// 1.0f / v;
    multData(dst);
  }
  
  public void checkMinMax(){
    themin = (float) 1e10;
    themax = (float) -1e10;
    for (int i = 0; i < numxyz; i++) {
      if (data[i] > themax) {
        themax = data[i];
      }
      if (data[i] < themin) {
        themin = data[i];
      }
    }  
  }
  
  public void datainvert() {

    float max = -1;
    for (int i = 0; i < numxyz; i++) {
      if(data[i]>max){
        max = data[i];
      }
    }
    
    for (int i = 0; i < numxyz; i++) {
      data[i] = max - data[i];
    }
    
  }
  
  public void datasubs(float v) {

    for (int i = 0; i < numxyz; i++) {
      data[i] = v - data[i];
    }
    
//    normalizeDataTo(1);
//    for (int i = 0; i < numxyz; i++) {
//      if(data[i] > 0)
//        data[i] = 0.0f;
//      else
//        data[i] = 1.0f;
//    }
  }

  public void post(){
    //System.out.print("\n");
    //for (int i = 0; i < numxyz; i++) {
    //  System.out.print(" "+data[i]);// = v - data[i];
    //}
    //System.out.print("\n");
  }
  
  
  
  // ///
  public void polygoniseData() {

    // Polygonise the grid
    // println("Polygonising data ...\n");
//    final int idx = 0;
    ntri = 0;
    trilist.clear();
    for (int i = 0; i < gx - 1; i++) {
      // if (i % (gx/10) == 0)
      // println("   Slice "+i+" of "+gx);
      for (int j = 0; j < gy - 1; j++) {
        for (int k = 0; k < gz - 1; k++) {
          grid.p[0].x = i * datastride.x;
          grid.p[0].y = j * datastride.y;
          grid.p[0].z = k * datastride.z;
          grid.val[0] = data[i + j * gx + k * gxgy];
          grid.p[1].x = (i + 1) * datastride.x;
          grid.p[1].y = j * datastride.y;
          grid.p[1].z = k * datastride.z;
          grid.val[1] = data[i + 1 + j * gx + k * gxgy];
          grid.p[2].x = (i + 1) * datastride.x;
          grid.p[2].y = (j + 1) * datastride.y;
          grid.p[2].z = k * datastride.z;
          grid.val[2] = data[i + 1 + (j + 1) * gx + k * gxgy];
          grid.p[3].x = i * datastride.x;
          grid.p[3].y = (j + 1) * datastride.y;
          grid.p[3].z = k * datastride.z;
          grid.val[3] = data[i + (j + 1) * gx + k * gxgy];
          grid.p[4].x = i * datastride.x;
          grid.p[4].y = j * datastride.y;
          grid.p[4].z = (k + 1) * datastride.z;
          grid.val[4] = data[i + j * gx + (k + 1) * gxgy];
          grid.p[5].x = (i + 1) * datastride.x;
          grid.p[5].y = j * datastride.y;
          grid.p[5].z = (k + 1) * datastride.z;
          grid.val[5] = data[i + 1 + j * gx + (k + 1) * gxgy];
          grid.p[6].x = (i + 1) * datastride.x;
          grid.p[6].y = (j + 1) * datastride.y;
          grid.p[6].z = (k + 1) * datastride.z;
          grid.val[6] = data[i + 1 + (j + 1) * gx + (k + 1) * gxgy];
          grid.p[7].x = i * datastride.x;
          grid.p[7].y = (j + 1) * datastride.y;
          grid.p[7].z = (k + 1) * datastride.z;
          grid.val[7] = data[i + (j + 1) * gx + (k + 1) * gxgy];
          final int n = Polygonise(grid, isolevel, triangles);

          // calc tri norms
          for (int a0 = 0; a0 < n; a0++) {
            triangles[a0].calcnormal(invertnormals);
          }


          for (int l = 0; l < n; l++) {
            final TRIANGLE t = new TRIANGLE(triangles[l]);
            trilist.add(t);
          }
          ntri += n;
        }
      }
    }
  }

  public void draw() {// PApplet applet) {

    for (int i = 0; i < trilist.size(); i++) {
      final TRIANGLE tri = trilist.get(i);
      applet.beginShape(PConstants.TRIANGLES);
      applet.normal((float)tri.n.x, (float)tri.n.y, (float)tri.n.z);
      applet.vertex(tri.p[0].x, tri.p[0].y, tri.p[0].z);
      applet.vertex(tri.p[1].x, tri.p[1].y, tri.p[1].z);
      applet.vertex(tri.p[2].x, tri.p[2].y, tri.p[2].z);
      applet.endShape();

      // triangle((float)tri[i].p[0].x, tri[i].p[0].y, tri[i].p[0].z,
      // tri[i].p[1].x, tri[i].p[1].y, tri[i].p[1].z,
      // tri[i].p[2].x, tri[i].p[2].y, tri[i].p[2].z);
    }
  }

  public void drawnormals(float s) {// PApplet applet, float s) {

    // println("draw ntri "+trilist.size());
    for (int i = 0; i < trilist.size(); i++) {
      final TRIANGLE tri = trilist.get(i);
      final float x = (tri.p[0].x + tri.p[1].x + tri.p[2].x) / 3.0f;
      final float y = (tri.p[0].y + tri.p[1].y + tri.p[2].y) / 3.0f;
      final float z = (tri.p[0].z + tri.p[1].z + tri.p[2].z) / 3.0f;

      applet.beginShape(PConstants.LINES);
      // normal((float)tri.n.x, (float)tri.n.y, (float)tri.n.z);
      applet.vertex(x, y, z);
      applet.vertex(x + tri.n.x * s, y + tri.n.y * s, z + tri.n.z * s);
      applet.endShape();

    }
  }

  
  public void drawsize(float s) {// PApplet applet) {

    
    float c = s / worldstride.x;///applet.map(s, istart, istop, ostart, ostop)
    
    // println("draw ntri "+trilist.size());
    for (int i = 0; i < trilist.size(); i++) {
      final TRIANGLE tri = trilist.get(i);
      applet.beginShape(PConstants.TRIANGLES);
      // normal((float)tri.n.x, (float)tri.n.y, (float)tri.n.z);
      applet.vertex(tri.p[0].x*c, tri.p[0].y*c, tri.p[0].z*c);
      applet.vertex(tri.p[1].x*c, tri.p[1].y*c, tri.p[1].z*c);
      applet.vertex(tri.p[2].x*c, tri.p[2].y*c, tri.p[2].z*c);
      applet.endShape();

    }
  }

  public void drawnormalssize(float siz, float s) {// PApplet applet, float s) {

    for (int i = 0; i < trilist.size(); i++) {
      final TRIANGLE tri = trilist.get(i);
      final float x = (tri.p[0].x + tri.p[1].x + tri.p[2].x) / 3.0f;
      final float y = (tri.p[0].y + tri.p[1].y + tri.p[2].y) / 3.0f;
      final float z = (tri.p[0].z + tri.p[1].z + tri.p[2].z) / 3.0f;

      applet.beginShape(PConstants.LINES);
      // normal((float)tri.n.x, (float)tri.n.y, (float)tri.n.z);
      applet.vertex(x, y, z);
      applet.vertex(x + tri.n.x * s, y + tri.n.y * s, z + tri.n.z * s);
      applet.endShape();

    }
  }
  
  
  public void drawGL(){
    
  }
  
  
  
  public void drawnormalsGL(float s){
    
  }
  

  /*
   * Given a grid cell and an isolevel, calculate the triangular facets
   * required to represent the isosurface through the cell. Return the number
   * of triangular facets, the array "triangles" will be loaded up with the
   * vertices at most 5 triangular facets. 0 will be returned if the grid cell
   * is either totally above of totally below the isolevel.
   */
  private int Polygonise(GRIDCELL grid, float isolevel, TRIANGLE triangles[]) {
    int i, ntriang;
    int cubeindex;
    final XYZ vertlist[] = new XYZ[12];
    /*
     * Determine the index into the edge table which tells us which vertices
     * are inside of the surface
     */
    cubeindex = 0;
    if (grid.val[0] < isolevel) {
      cubeindex |= 1;
    }
    if (grid.val[1] < isolevel) {
      cubeindex |= 2;
    }
    if (grid.val[2] < isolevel) {
      cubeindex |= 4;
    }
    if (grid.val[3] < isolevel) {
      cubeindex |= 8;
    }
    if (grid.val[4] < isolevel) {
      cubeindex |= 16;
    }
    if (grid.val[5] < isolevel) {
      cubeindex |= 32;
    }
    if (grid.val[6] < isolevel) {
      cubeindex |= 64;
    }
    if (grid.val[7] < isolevel) {
      cubeindex |= 128;
    }

    /* Cube is entirely in/out of the surface */
    if (edgeTable[cubeindex] == 0) {
      return (0);
    }

    /* Find the vertices where the surface intersects the cube */
    // int temp = edgeTable[cubeindex] & 1;
    if ((edgeTable[cubeindex] & 1) != 0) {
      vertlist[0] = VertexInterp(isolevel, grid.p[0], grid.p[1],
          grid.val[0], grid.val[1]);
    }
    if ((edgeTable[cubeindex] & 2) != 0) {
      vertlist[1] = VertexInterp(isolevel, grid.p[1], grid.p[2],
          grid.val[1], grid.val[2]);
    }
    if ((edgeTable[cubeindex] & 4) != 0) {
      vertlist[2] = VertexInterp(isolevel, grid.p[2], grid.p[3],
          grid.val[2], grid.val[3]);
    }
    if ((edgeTable[cubeindex] & 8) != 0) {
      vertlist[3] = VertexInterp(isolevel, grid.p[3], grid.p[0],
          grid.val[3], grid.val[0]);
    }
    if ((edgeTable[cubeindex] & 16) != 0) {
      vertlist[4] = VertexInterp(isolevel, grid.p[4], grid.p[5],
          grid.val[4], grid.val[5]);
    }
    if ((edgeTable[cubeindex] & 32) != 0) {
      vertlist[5] = VertexInterp(isolevel, grid.p[5], grid.p[6],
          grid.val[5], grid.val[6]);
    }
    if ((edgeTable[cubeindex] & 64) != 0) {
      vertlist[6] = VertexInterp(isolevel, grid.p[6], grid.p[7],
          grid.val[6], grid.val[7]);
    }
    if ((edgeTable[cubeindex] & 128) != 0) {
      vertlist[7] = VertexInterp(isolevel, grid.p[7], grid.p[4],
          grid.val[7], grid.val[4]);
    }
    if ((edgeTable[cubeindex] & 256) != 0) {
      vertlist[8] = VertexInterp(isolevel, grid.p[0], grid.p[4],
          grid.val[0], grid.val[4]);
    }
    if ((edgeTable[cubeindex] & 512) != 0) {
      vertlist[9] = VertexInterp(isolevel, grid.p[1], grid.p[5],
          grid.val[1], grid.val[5]);
    }
    if ((edgeTable[cubeindex] & 1024) != 0) {
      vertlist[10] = VertexInterp(isolevel, grid.p[2], grid.p[6],
          grid.val[2], grid.val[6]);
    }
    if ((edgeTable[cubeindex] & 2048) != 0) {
      vertlist[11] = VertexInterp(isolevel, grid.p[3], grid.p[7],
          grid.val[3], grid.val[7]);
    }

    /* Create the triangle */
    ntriang = 0;
    for (i = 0; triTable[cubeindex][i] != -1; i += 3) {
      triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i]];
      triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i + 1]];
      triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i + 2]];
      ntriang++;
    }

    return (ntriang);
  }

  /*
   * Linearly interpolate the position where an isosurface cuts an edge
   * between two vertices, each with their own scalar value
   */
  private XYZ VertexInterp(float isolevel, XYZ p1, XYZ p2, float valp1,
      float valp2) {
    float mu;
    final XYZ p = new XYZ();
    p.x = p.y = p.z = 0.0f;

    if (Math.abs(isolevel - valp1) < 0.00001) {
      return (p1);
    }
    if (Math.abs(isolevel - valp2) < 0.00001) {
      return (p2);
    }
    if (Math.abs(valp1 - valp2) < 0.00001) {
      return (p1);
    }
    mu = ((isolevel - valp1) / (valp2 - valp1));
    p.x = p1.x + mu * (p2.x - p1.x);
    p.y = p1.y + mu * (p2.y - p1.y);
    p.z = p1.z + mu * (p2.z - p1.z);

    return (p);
  }

  // /TABLES
  int edgeTable[] = { // 256
  0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f,
      0xb06, 0xc0a, 0xd03, 0xe09, 0xf00, 0x190, 0x99, 0x393, 0x29a,
      0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895, 0xb9f, 0xa96, 0xd9a,
      0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x33, 0x13a, 0x636, 0x73f,
      0x435, 0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39,
      0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa, 0x7a6, 0x6af, 0x5a5, 0x4ac,
      0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 0x460,
      0x569, 0x663, 0x76a, 0x66, 0x16f, 0x265, 0x36c, 0xc6c, 0xd65,
      0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3,
      0x6fa, 0x1f6, 0xff, 0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff, 0xef6,
      0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a, 0x256,
      0x35f, 0x55, 0x15c, 0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53,
      0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5,
      0xcc, 0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
      0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc,
      0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0, 0x950, 0x859,
      0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55, 0x35f,
      0x256, 0x55a, 0x453, 0x759, 0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa,
      0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc, 0x3f5, 0xff, 0x1f6, 0x6fa,
      0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f,
      0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x66, 0x76a, 0x663, 0x569,
      0x460, 0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
      0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa, 0x1a3, 0x2a9, 0x3a0, 0xd30,
      0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435,
      0x73f, 0x636, 0x13a, 0x33, 0x339, 0x230, 0xe90, 0xf99, 0xc93,
      0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795, 0x49f, 0x596,
      0x29a, 0x393, 0x99, 0x190, 0xf00, 0xe09, 0xd03, 0xc0a, 0xb06,
      0xa0f, 0x905, 0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203,
      0x109, 0x0 };
  int triTable[][] = // 256x16
  { { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
      { 2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
      { 8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1 },
      { 3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1 },
      { 4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
      { 4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
      { 5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
      { 2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1 },
      { 9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
      { 2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1 },
      { 10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1 },
      { 5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
      { 5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1 },
      { 10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1 },
      { 8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1 },
      { 2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
      { 7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
      { 2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1 },
      { 11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1 },
      { 5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1 },
      { 11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1 },
      { 11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
      { 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1 },
      { 2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
      { 5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1 },
      { 6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
      { 3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1 },
      { 6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
      { 5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
      { 10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
      { 6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1 },
      { 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1 },
      { 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1 },
      { 3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
      { 5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
      { 0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1 },
      { 9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1 },
      { 8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
      { 5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1 },
      { 0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1 },
      { 6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1 },
      { 10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1 },
      { 10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
      { 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
      { 1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1 },
      { 0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
      { 10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1 },
      { 3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
      { 6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1 },
      { 9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1 },
      { 8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1 },
      { 3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
      { 6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1 },
      { 10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1 },
      { 10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
      { 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1 },
      { 7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
      { 7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
      { 2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1 },
      { 1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1 },
      { 11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1 },
      { 8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1 },
      { 0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1 },
      { 7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
      { 10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1 },
      { 2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1 },
      { 6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1 },
      { 7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
      { 2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1 },
      { 10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
      { 10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1 },
      { 0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1 },
      { 7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
      { 6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
      { 8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1 },
      { 6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1 },
      { 4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1 },
      { 10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1 },
      { 8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1 },
      { 1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
      { 8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1 },
      { 10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1 },
      { 10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
      { 5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },
      { 11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1 },
      { 9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },
      { 6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1 },
      { 7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1 },
      { 3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1 },
      { 7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1 },
      { 3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1 },
      { 6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1 },
      { 9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1 },
      { 1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1 },
      { 4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1 },
      { 7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1 },
      { 6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1 },
      { 0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1 },
      { 6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1 },
      { 0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1 },
      { 11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1 },
      { 6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1 },
      { 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1 },
      { 9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1 },
      { 1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1 },
      { 10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1 },
      { 0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1 },
      { 5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1 },
      { 10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1 },
      { 11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1 },
      { 9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1 },
      { 7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1 },
      { 2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
      { 8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1 },
      { 9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1 },
      { 9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1 },
      { 1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1 },
      { 5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1 },
      { 0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1 },
      { 10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1 },
      { 2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1 },
      { 0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1 },
      { 0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1 },
      { 9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1 },
      { 5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1 },
      { 5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1 },
      { 8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1 },
      { 9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1 },
      { 1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1 },
      { 3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1 },
      { 4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1 },
      { 9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1 },
      { 11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
      { 11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1 },
      { 2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1 },
      { 9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1 },
      { 3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1 },
      { 1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1 },
      { 4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1 },
      { 0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
      { 9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1 },
      { 1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
      { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 } };

  /*
   * 
   * 
   * file io
   */

  
  
  
  /**
   * return the version of the library.
   * 
   * @return String
   */
  public String version() {
  
    print("Version: 10");
    return "10";

  }
  
}


public class TRIANGLE {
  public XYZ p[] = new XYZ[3];
  public XYZ n = new XYZ(0, 1, 0);

  public TRIANGLE() {
    for (int i = 0; i < 3; i++) {
      p[i] = new XYZ();
    }
  }

  public TRIANGLE(TRIANGLE t) {
    n = new XYZ(t.n);
    for (int i = 0; i < 3; i++) {
      p[i] = new XYZ(t.p[i]);
    }
  }

  public void calcnormal(boolean invertnormals) {
    // PVector nv = new PVector();
    // PVector v2v1 = new PVector( (float)(p[1].x-p[0].x), (float)(p[1].y -
    // p[0].y), (float)(p[1].z - p[0].z ) );
    // PVector v3v1 = new PVector( (float)(p[2].x-p[0].x), (float)(p[2].y -
    // p[0].y), (float)(p[2].z - p[0].z) );

    XYZ nv = new XYZ();
    XYZ v1 = new XYZ((p[1].x - p[0].x), (p[1].y - p[0].y),
        (p[1].z - p[0].z));
    XYZ v2 = new XYZ((p[2].x - p[0].x), (p[2].y - p[0].y),
        (p[2].z - p[0].z));

    nv.x = (v1.y * v2.z) - (v1.z * v2.y);
    nv.y = (v1.z * v2.x) - (v1.x * v2.z);
    nv.z = (v1.x * v2.y) - (v1.y * v2.x);
    // recheck, this
    // nv.normalize();
    // float d = nv.x*nv.x+nv.y*nv.y+nv.z*nv.z;
    // if (d>0) {
    // d = (float) (1.0/Math.sqrt(d));
    // nv.x*=d;
    // nv.y*=d;
    // nv.z*=d;
    // }
    if (!invertnormals) {
      n.x = nv.x;
      n.y = nv.y;
      n.z = nv.z;
    } else {
      n.x = -nv.x;
      n.y = -nv.y;
      n.z = -nv.z;
    }
  }
}



public class XYZ {
 public float x, y, z;
 public XYZ() {
  }
 public XYZ(float a, float b, float c) {
    x=a; 
    y=b; 
    z=c;
  }
 public XYZ(XYZ a) {
    x=a.x; 
    y=a.y; 
    z=a.z;
  }
}




class Cell {
  private boolean activated;
  private boolean hidden;
  private boolean aboveCell;
  private float concentration;
  public float r1, r2, c;
  public Cell(float r1, float r2, float c, boolean active) {
    activated = active;
    hidden = false;
    aboveCell = false;
    concentration = 0;
    this.r1 = r1;
    this.r2 = r2;
    this.c = c;
  }
  public Cell(Cell c) {
    activated = c.activated;
    hidden = false;
    aboveCell = false;
    concentration = c.concentration;
    r1 = c.r1;
    r2 = c.r2;
    this.c = c.c;
  }
  public void activate() {
    activated = true;
  }
  public void deactivate() {
    activated = false;
  }
  public boolean isActivated() {
    return activated;
  }
  public void addConcentration(float c) {
    concentration += c;
  }
  public void step() {
    activated = concentration > 0;
    concentration = 0;
  }
  public boolean isHidden() {
    return hidden;
  }
  public boolean aboveCell() {
    return aboveCell;
  }
  public void markAboveCell() {
    aboveCell = true;
  }
  public void markHidden() {
    hidden = true;
  }
  public void setR1(float val) {
    r1 = val;
  }
  public void setR2(float val) {
    r2 = val;
  }
}

class Grid {
  Cell[][] grid;
  Cell[][] newGrid;
  int w, h;
  float r1g, r2g, c;
  public Grid(int w, int h, float r1, float r2, float c) {
    this.w = w;
    this.h = h;
    this.r1g = r1;
    this.r2g = r2;
    this.c = c;
    grid = new Cell[w][h];
    newGrid = new Cell[w][h];
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        //grid[i][j] = new Cell(lerp(r1-2, r1+2, (float)i/w), lerp(r2-2, r2+2, (float)j/h), c);
        grid[i][j] = new Cell(r1, r2, c, (abs(float(i)/w-.5)<.1 && abs(float(j)/h-.5)<.1 && random(1)<.5));
        newGrid[i][j] = new Cell(grid[i][j]);
      }
    }
  }
  public void distribute(Cell[][] target, float conc, int x, int y, float r) {
    for (int i0=int(x-r); i0<=x+r; i0++) {
      for (int j0=int(y-r); j0<=y+r; j0++) {
        float d = sqrt((x-i0)*(x-i0)+(y-j0)*(y-j0));
        if (d <= r) {
          if (i0 < 0 || i0 >= w || j0 < 0 || j0 >= h) {
          }
          else {
            target[i0][j0].addConcentration(d/r*conc);
          }
        }
      }
    }
  }
  public Cell getCell(int i, int j) {
    if (i >= 0 && i < w && j >= 0 && j < h) {
      return grid[i][j];
    }
    else {
      return new Cell(r1g, r2g, c, false);
    }
  }
  public boolean activatedOnGrid(int i, int j) {
    if (i >= 0 && i < w && j >= 0 && j < h) {
      return grid[i][j].isActivated();
    }
    else {
      return false;
    }
  }
  public void updateNewGrid() {
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        if (activatedOnGrid(i, j)) {
          distribute(newGrid, 1, i, j, grid[i][j].r1);
          distribute(newGrid, grid[i][j].c, i, j, grid[i][j].r2);
        }
      }
    }
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        newGrid[i][j].step();
      }
    }
  }
  public void markHidden() {
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        if (newGrid[i][j].isActivated() && grid[i][j].isActivated()) {
          newGrid[i][j].markAboveCell();
          if (activatedOnGrid(i-1, j) && activatedOnGrid(i+1, j) && activatedOnGrid(i, j-1) && activatedOnGrid(i, j+1)) {          
            grid[i][j].markHidden();
          }
        }
      }
    }
  }
  public void switchToNewGrid() {
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        grid[i][j] = new Cell(newGrid[i][j]);
      }
    }
  }
  
  public void draw(float thickness) {
    noFill();
    for (float i=0; i<w; i++) {
      for (float j=0; j<h; j++) {
        if (activatedOnGrid(int(i), int(j))) {
          stroke(0);
        }
        else {
          noStroke();
        }
        pushMatrix();
        translate(i/w*width, j/h*height);
        box(width/w, height/h, thickness);
        popMatrix();
      }
    }
  }
  public void setNewRadii(float r01, float r02) {
    for (int i=0; i<w; i++) {
      for (int j=0; j<h; j++) {
        Cell cel = grid[i][j];
        cel.setR1(r01);
        cel.setR2(r02);
      }
    }
  }
}

class GridBuilder {
  private Grid g;
  float z;
  float dz;
  ArrayList<GridLayer> layers;
  
  class GridLayer {
    ArrayList<PVector> boxes;
    int w, h;
    public GridLayer(Grid mygrid, float z) {
      boxes = new ArrayList<PVector>();
      w = mygrid.w;
      h = mygrid.h;
      for (int i=0; i<g.w; i++) {
        for (int j=0; j<g.h; j++) {
          if (mygrid.activatedOnGrid(i,j) && !mygrid.getCell(i, j).isHidden()) {
            boxes.add(new PVector((float)i/mygrid.w*width, (float)j/mygrid.h*height, z));
          }
        }
      }
    }
    public void draw(float thickness) {
      for (PVector b : boxes) {
        pushMatrix();
        translate(b.x, b.y, b.z);
        box(width/w, height/h, thickness);
        popMatrix();
      }
    }
    public int countBoxes() {
      return boxes.size();
    }
  }
  
  public GridBuilder(float z0, float dz, int w, int h, float r11, float r22, float c) {
    g = new Grid(w, h, r11, r22, c);
    this.z = z0;
    this.dz = dz;
    layers = new ArrayList<GridLayer>();
  }
  public void update() {
    // println("Starting update");
    g.updateNewGrid();
    // println("Updated new grid");
    if (layers.size() > 0) {
      g.markHidden();
      layers.set(layers.size()-1, new GridLayer(g, z - dz));
    }
    // println("Hid hidden blocks");
    g.switchToNewGrid();
    // println("Switched to new grid");
    layers.add(new GridLayer(g, z));
    // println("Added new layer to layers");
    z += dz;
    // println("Incremented z");
  }
  public void updateToMarchingCubes(MarchingCubes mcubes) {
    g.updateNewGrid();
    g.markHidden();
    g.switchToNewGrid();
    for (int i=0; i<g.w; i++) {
      for (int j=0; j<g.h; j++) {
        if (g.activatedOnGrid(i,j) && !g.getCell(i, j).isHidden()) {
          mcubes.addIsoPoint(isoval, new XYZ((float)i/g.w * width, (float)j/g.h * height, z));
        }
      }
    }
    z += dz;
  }
  public void draw(float thickness) {
    for (GridLayer l : layers) {
      l.draw(thickness);
    }
  }
  
  public Grid getGrid() {
    return g;
  }
  public void reset(float z0, float are1, float are2) {
    z = z0;
    g = new Grid(g.w, g.h, are1, are2, g.c);
    layers = new ArrayList<GridLayer>();
  }
}


MarchingCubes mc;
GridBuilder gb;
String folder;
float activatorRad = 4;
float inhibitorRad = 8;
float c = -.25;

PVector rot;

float isoval = .02;
int mcres = 130;

void setup() {
  size(900, 900, P3D);
  //fullScreen(P3D);
  folder = month() + "-" + day() + "-" + year() + "/" + hour() + ";" + minute() + "," + second();
  gb = new GridBuilder(-500, 5, 200, 200, activatorRad, inhibitorRad, c);
  //gb.update();
  background(255);
  
  mc = new MarchingCubes(this, width, height, -500, mcres, mcres, mcres);    
  mc.isolevel = 0.0161;
  mc.zeroData();
  
  mouseClicked();
  
  rot = new PVector(0, 0);
}

void draw() {
  
  noStroke();
  fill(128);
  
  //lights();
  
  pushMatrix();
  
  //background(240);
  
  translate(width/2, height/2, gb.z);
  
  rotateY(rot.x);
  rotateX(rot.y);
  
  translate(-width/2, -height/2, -gb.z);
  
  //gb.draw(5);
  
  lights();
  background(0);
  
  mc.draw();
  //mc.drawnormals(0.01);
  
  popMatrix();
  
  noLights();
  fill(255);
  text("Press left-click or hold right-click to simulate.\nSpace to restart", 10, 10);
  
  if (mousePressed && mouseButton == RIGHT) {
    if (frameCount % 2 == 0) {
      mouseClicked();
    }
  }
}

void mouseClicked() {
  gb.updateToMarchingCubes(mc);
  mc.polygoniseData();
  
}

void mouseDragged() {
  if (mouseButton==LEFT) {
      rot.x += -(mouseX-pmouseX)*0.01;
      rot.y += (mouseY-pmouseY)*0.01;
    }
}

void keyPressed() {
  if (key == ' ') {
    // println("Resetting");
    //gb = new GridBuilder(-500, 5, 100, 100, activatorRad, inhibitorRad, c);
    mc.zeroData();
    gb.reset(-500, activatorRad, inhibitorRad);
    mouseClicked();
    // println("Finished re-initializing");
  }
}

void setActivatorRadius(float r) {
  activatorRad = r;
  gb.getGrid().setNewRadii(activatorRad, inhibitorRad);
  // println("Set new activator radius to " + activatorRad);
}

void setInhibitorRadius(float r) {
  inhibitorRad = r;
  gb.getGrid().setNewRadii(activatorRad, inhibitorRad);
  // println("Set new inhibitor radius to " + inhibitorRad);
}