BEGIN VERTEX
void main()
{
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position.xyz, 1.0);
}
END VERTEX

BEGIN FRAGMENT
const int NUM_LIGHTS = 4;
const int NUM_OBJECTS = 2;
//--------------------------------------------
//Global Constants
//--------------------------------------------
const int MAX_MARCHING_STEPS = 127;
const float MIN_DIST = 0.0;
const float EPSILON = 0.0001;
const vec4 ORIGIN = vec4(0,0,0,1);
//--------------------------------------------
//Global Variables
//--------------------------------------------
vec4 sampleEndPoint = vec4(1, 1, 1, 1);
vec4 sampleTangentVector = vec4(1, 1, 1, 1);
mat4 totalFixMatrix = mat4(1.0);
//normal vector
vec4 N = ORIGIN; 
vec4 globalLightColor = ORIGIN;
int hitWhich = 0;
//-------------------------------------------
//Translation & Utility Variables
//--------------------------------------------
uniform int isStereo;
//uniform int geometry;
uniform vec2 screenResolution;
uniform float fov;
uniform mat4 invGenerators[6];
uniform mat4 currentBoost;
uniform mat4 stereoBoosts[2];
uniform mat4 cellBoost; 
uniform mat4 invCellBoost;
uniform int maxSteps;
uniform float maxDist;
//--------------------------------------------
//Lighting Variables & Global Object Variables
//--------------------------------------------
uniform vec4 lightPositions[4];
//w component is the light's attenuation -- 6 since we need controllers
uniform vec4 lightIntensities[6]; 
uniform int attnModel;
uniform sampler2D texture;
//Max is two
uniform int controllerCount; 
uniform mat4 controllerBoosts[2];
//uniform vec4 controllerDualPoints[6];
uniform mat4 globalObjectBoosts[4];
uniform mat4 invGlobalObjectBoosts[4]; 
uniform vec3 globalObjectRadii[4];
uniform int globalObjectTypes[4];
//--------------------------------------------
//Scene Dependent Variables
//--------------------------------------------
uniform vec4 halfCubeDualPoints[3];
uniform float halfCubeWidthKlein;
uniform float sphereRad;
uniform float tubeRad;
uniform vec4 vertexPosition;
uniform float vertexSurfaceOffset;

// These are the planar mirrors of the fundamental simplex in the Klein (or analagous) model.
// Order is mirrors opposite: vertex, edge, face, cell.
// The xyz components of a vector give the unit normal of the mirror. The sense will be that the normal points to the outside of the simplex.
// The w component is the offset from the origin.
uniform bool useSimplex;
uniform vec4 simplexMirrorsKlein[4];

// The type of cut (1=sphere, 2=horosphere, 3=plane) for the vertex opposite the fundamental simplex's 4th mirror.
// These integers match our values for the geometry of the honeycomb vertex figure.
// We'll need more of these later when we support more symmetry groups.
uniform int cut4;

//Raymarch Functions
float unionSDF(float d1, float d2){
  return min(d1, d2);
}

//--------------------------------------------------------------------
// Hyperbolic Functions
//--------------------------------------------------------------------
float cosh(float x){
  float eX = exp(x);
  return (0.5 * (eX + 1.0/eX));
}

float acosh(float x){ //must be more than 1
  return log(x + sqrt(x*x-1.0));
}

//--------------------------------------------------------------------
// Euclidean Functions
//--------------------------------------------------------------------

float eucDot(vec4 u, vec4 v){
  return dot(u.xyz,v.xyz);
}

float eucNorm(vec4 v){
  return sqrt(abs(eucDot(v,v)));
}

vec4 eucNormalize(vec4 u, bool toTangent){
  if(toTangent){
    u.xyz = normalize(u.xyz);
    u.w = 0.0;
    return u;
  }
  else{
    u.w = 1.0;
    return u;
  }
}

float eucDistance(vec4 u, vec4 v){
  return distance(u.xyz, v.xyz);
}

//Given two positions find the unit tangent vector at the first that points to the second
vec4 eucDirection(vec4 u, vec4 v){ 
  vec4 w = v-u;
  return eucNormalize(w, true);
}

//calculate the new direction vector (v) for the continuation of the ray from the new ray origin (u)
//having moved by fix matrix
vec4 eucFixDirection(vec4 u, vec4 v, mat4 fixMatrix){
  return v;
}

vec4 projectToKlein(vec4 v)
{
  // We are already effectively Klein (i.e. lines are straight in the model)
  v.w = 1.0;
  return v;
}

vec4 pointOnGeodesic(vec4 u, vec4 vPrime, float dist)
{ 
  // get point at distance dist on the geodesic from u in the direction vPrime
  return projectToKlein( u + vPrime*dist );
}

vec4 tangentVectorOnGeodesic(vec4 u, vec4 vPrime, float dist)
{
  return vPrime;
}

vec4 pointOnGeodesicAtInfinity(vec4 u, vec4 vPrime)
{ 
  // I'm not yet sure what we should be doing here.
  return projectToKlein(u + vPrime);
}

//--------------------------------------------------------------------
// SDFs
//--------------------------------------------------------------------

float sphereSDF(vec4 samplePoint, vec4 center, float radius){
  return eucDistance(samplePoint, center) - radius;
}

float sortOfEllipsoidSDF(vec4 samplePoint, mat4 boostMatrix){
  return sphereSDF(eucNormalize(samplePoint * boostMatrix, false), ORIGIN, 0.05);
}

float controllerSDF(vec4 samplePoint, mat4 controllerBoost, float radius){
  float sphere = sphereSDF(samplePoint, ORIGIN * controllerBoost, radius);
  
  //generated on JS side
  //may be subject to change
  //translateByVector(0,0,0.2)*scale matrix (0.8, 0.8, 0.4)
  mat4 scaleMatrix = mat4(
    0.8, 0.0, 0.0, 0.0,
    0.0, 0.8, 0.0, 0.0,
    0.0, 0.0, 0.408, 0.0805,
    0.0, 0.0, 0.2013, 1.02
  );

  //We need to offset this so that the ellipsoid is not centered at the same point as the sphere
  float ellipsoid = sortOfEllipsoidSDF(samplePoint, scaleMatrix * controllerBoost);
  return unionSDF(sphere, ellipsoid);
}

float geodesicCylinderHSDFplanes(vec4 samplePoint, vec4 cylinderCorePoint, vec4 direction, float radius)
{
  cylinderCorePoint = halfCubeDualPoints[0] + halfCubeDualPoints[1];
  direction = halfCubeDualPoints[2];

  vec4 pos = (samplePoint - cylinderCorePoint);
  return length(pos.xyz - eucDot(pos, direction) * direction.xyz) - radius;
}

//
// Functions below this don't apply, but we need them included to make the shader compile.
//
float horosphereHSDF(vec4 samplePoint, vec4 lightPoint, float offset)
{
  return 0.0;
}

float geodesicPlaneHSDF(vec4 samplePoint, vec4 dualPoint, float offset){
  return 0.0;
}

//--------------------------------------------------------------------
// Scene Functions
//--------------------------------------------------------------------

float localSceneSDF(vec4 samplePoint){
    float sphere = sphereSDF(samplePoint, ORIGIN, sphereRad);
    float vertexSphere = 0.0;
    if(cut4 == 1) {
        vertexSphere = sphereSDF(abs(samplePoint), vertexPosition, vertexSurfaceOffset);
    }
    else if(cut4 == 2) {
        vertexSphere = horosphereHSDF(abs(samplePoint), vertexPosition, vertexSurfaceOffset);
    }
    else if(cut4 == 3) {
        vertexSphere = geodesicPlaneHSDF(abs(samplePoint), vertexPosition, vertexSurfaceOffset);
    }
    float final = -unionSDF(vertexSphere,sphere);
    return final;
}

//GLOBAL OBJECTS SCENE ++++++++++++++++++++++++++++++++++++++++++++++++
float globalSceneSDF(vec4 samplePoint){
  vec4 absoluteSamplePoint = samplePoint * cellBoost; // correct for the fact that we have been moving
  float distance = maxDist;
  //Light Objects
  for(int i=0; i<NUM_LIGHTS; i++){
    float objDist;
    if(lightIntensities[i].w == 0.0) { objDist = maxDist; }
    else{
      objDist = sphereSDF(absoluteSamplePoint, lightPositions[i], 1.0/(10.0*lightIntensities[i].w));
      distance = min(distance, objDist);
      if(distance < EPSILON){
        hitWhich = 1;
        globalLightColor = lightIntensities[i];
        return distance;
      }
    }
  }
  //Controller Objects
  for(int i=0; i<2; i++){
    if(controllerCount != 0){
      float objDist = sphereSDF(samplePoint, ORIGIN*controllerBoosts[i]*currentBoost, 1.0/(10.0 * lightIntensities[i+NUM_LIGHTS].w));
      //float objDist = controllerSDF(absoluteSamplePoint, controllerBoosts[i]*currentBoost, 1.0/(10.0 * lightIntensities[i+NUM_LIGHTS].w));
      distance = min(distance, objDist);
      if(distance < EPSILON){
        hitWhich = 1;
        globalLightColor = lightIntensities[i+4];
        return distance;
      }
      if(controllerCount == 1) break;
    }
  }
  //Global Objects
  for(int i=0; i<NUM_OBJECTS; i++) {
    float objDist;
    if(length(globalObjectRadii[i]) == 0.0){ objDist = maxDist;}
    else{
      if(globalObjectTypes[i] == 0) { objDist = sphereSDF(eucNormalize(absoluteSamplePoint * globalObjectBoosts[i], false), ORIGIN, globalObjectRadii[i].x); }
      else if(globalObjectTypes[i] == 1) { objDist = sortOfEllipsoidSDF(absoluteSamplePoint, globalObjectBoosts[i]);}
      else { objDist = maxDist; }
      distance = min(distance, objDist);
      if(distance < EPSILON){
        hitWhich = 2;
      }
    }
  }
  return distance;
}
  
//--------------------------------------------------------------------
// Lighting Functions
//--------------------------------------------------------------------
vec4 texcube(vec4 samplePoint, mat4 toOrigin){
    float k = 4.0;
    vec4 newSP = samplePoint * toOrigin;
    vec3 p = mod(newSP.xyz,1.0);
    vec3 n = eucNormalize(N*toOrigin, true).xyz; //Very hacky you are warned
    vec3 m = pow(abs(n), vec3(k));
    vec4 x = texture2D(texture, p.yz);
    vec4 y = texture2D(texture, p.zx);
    vec4 z = texture2D(texture, p.xy);
    return (x*m.x + y*m.y + z*m.z) / (m.x+m.y+m.z);
}

float attenuation(float distToLight, vec4 lightIntensity){
  float att;
  if(attnModel == 1) //Inverse Linear
    att  = 0.75/ (0.01+lightIntensity.w * distToLight);  
  else if(attnModel == 2) //Inverse Square
    att  = 1.0/ (0.01+lightIntensity.w * distToLight* distToLight);
  else if(attnModel == 3) // Inverse Cube
    att = 1.0/ (0.01+lightIntensity.w*distToLight*distToLight*distToLight);
  else if(attnModel == 4) //Physical
    att  = 1.0/ (0.01+lightIntensity.w*cosh(2.0*distToLight)-1.0);
  else //None
    att  = 0.25; //if its actually 1 everything gets washed out
  return att;
}

vec3 lightingCalculations(vec4 SP, vec4 TLP, vec4 V, vec3 baseColor, vec4 lightIntensity){
  float distToLight = eucDistance(SP, TLP);
  float att = attenuation(distToLight, lightIntensity);
  //Calculations - Phong Reflection Model
  float shadow = 1.0;
  vec4 L = eucDirection(SP, TLP);
  vec4 R = 2.0*eucDot(L, N)*N - L;
  //Calculate Diffuse Component
  float nDotL = max(eucDot(N, L),0.0);
  vec3 diffuse = lightIntensity.rgb * nDotL;
  //Calculate Specular Component
  float rDotV = max(eucDot(R, V),0.0);
  vec3 specular = lightIntensity.rgb * pow(rDotV,10.0);
  //Compute final color
  return att*(shadow*((diffuse*baseColor) + specular));
}

vec3 phongModel(mat4 invObjectBoost, bool isGlobal){
    vec4 V, samplePoint;
    float ambient = 0.1;
    vec3 baseColor = vec3(0.0,1.0,1.0);

    //--------------------------------------------
    //Setup Variables
    //--------------------------------------------    
	samplePoint = sampleEndPoint;
    V = -sampleTangentVector; //Viewer is in the direction of the negative ray tangent vector
    if(isGlobal){ //this may be possible to move outside function as we already have an if statement for global v. local
      totalFixMatrix = mat4(1.0);
      baseColor = texcube(samplePoint, cellBoost * invObjectBoost).xyz; 
    }
    else{
      baseColor = texcube(samplePoint, mat4(1.0)).xyz;
    }

    //Setup up color with ambient component
    vec3 color = baseColor * ambient; 

    //--------------------------------------------
    //Lighting Calculations
    //--------------------------------------------
    vec4 translatedLightPosition;
    //Standard Light Objects
    for(int i = 0; i<NUM_LIGHTS; i++){
      if(lightIntensities[i].w != 0.0){
        translatedLightPosition = lightPositions[i]*invCellBoost*totalFixMatrix;
        color += lightingCalculations(samplePoint, translatedLightPosition, V, baseColor, lightIntensities[i]);
      }
    }
    //Lights for Controllers
    for(int i = 0; i<2; i++){
      if(controllerCount == 0) break; //if there are no controllers do nothing
      else translatedLightPosition = ORIGIN*controllerBoosts[i]*currentBoost*totalFixMatrix;
      color += lightingCalculations(samplePoint, translatedLightPosition, V, baseColor, lightIntensities[i+4]);
      if(controllerCount == 1) break; //if there is one controller only do one loop
    }
    return color;
}

//NORMAL FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++
vec4 estimateNormal(vec4 p) { // normal vector is in tangent hyperplane to hyperboloid at p
    // float denom = sqrt(1.0 + p.x*p.x + p.y*p.y + p.z*p.z);  // first, find basis for that tangent hyperplane
    float newEp = EPSILON * 10.0;
    vec4 basis_x = eucNormalize(vec4(p.w,0.0,0.0,p.x), true);  // dw/dx = x/w on hyperboloid
    vec4 basis_y = vec4(0.0,p.w,0.0,p.y);  // dw/dy = y/denom
    vec4 basis_z = vec4(0.0,0.0,p.w,p.z);  // dw/dz = z/denom  /// note that these are not orthonormal!
    basis_y = eucNormalize(basis_y - eucDot(basis_y, basis_x)*basis_x, true); // need to Gram Schmidt
    basis_z = eucNormalize(basis_z - eucDot(basis_z, basis_x)*basis_x - eucDot(basis_z, basis_y)*basis_y, true);
    if(hitWhich == 1 || hitWhich == 2){ //global light scene
      return eucNormalize( //p+EPSILON*basis_x should be lorentz normalized however it is close enough to be good enough
          basis_x * (globalSceneSDF(p + newEp*basis_x) - globalSceneSDF(p - newEp*basis_x)) +
          basis_y * (globalSceneSDF(p + newEp*basis_y) - globalSceneSDF(p - newEp*basis_y)) +
          basis_z * (globalSceneSDF(p + newEp*basis_z) - globalSceneSDF(p - newEp*basis_z)),
          true
      );
    }
    else{ //local scene
      return eucNormalize(
          basis_x * (localSceneSDF(p + newEp*basis_x) - localSceneSDF(p - newEp*basis_x)) +
          basis_y * (localSceneSDF(p + newEp*basis_y) - localSceneSDF(p - newEp*basis_y)) +
          basis_z * (localSceneSDF(p + newEp*basis_z) - localSceneSDF(p - newEp*basis_z)),
          true
      );
    }
}

//--------------------------------------------------------------------
// Marching Functions
//--------------------------------------------------------------------

vec4 getRayPoint(vec2 resolution, vec2 fragCoord, bool isLeft){ //creates a point that our ray will go through
    if(isStereo == 1){
      resolution.x = resolution.x * 0.5;
      if(!isLeft) { fragCoord.x = fragCoord.x - resolution.x; }
    }
    vec2 xy = 0.2*((fragCoord - 0.5*resolution)/resolution.x);
    float z = 0.1/tan(radians(fov*0.5));
    vec4 p =  eucNormalize(vec4(xy,-z,1.0), false);
    return p;
}

bool isOutsideSimplex(vec4 samplePoint, out mat4 fixMatrix){
  vec4 kleinSamplePoint = projectToKlein(samplePoint);
  for(int i=0; i<4; i++){
    vec3 normal = simplexMirrorsKlein[i].xyz;
    vec3 offsetSample = kleinSamplePoint.xyz - normal * simplexMirrorsKlein[i].w;  // Deal with any offset.
    if( dot(offsetSample, normal) > 1e-7 ) {
      fixMatrix = invGenerators[i];
      return true;
    }
  }
  return false;
}

// This function is intended to be geometry-agnostic.
bool isOutsideCell(vec4 samplePoint, out mat4 fixMatrix){
  if( useSimplex ) {
    return isOutsideSimplex( samplePoint, fixMatrix );
  }

  vec4 kleinSamplePoint = projectToKlein(samplePoint);
  if(kleinSamplePoint.x > halfCubeWidthKlein){
    fixMatrix = invGenerators[0];
    return true;
  }
  if(kleinSamplePoint.x < -halfCubeWidthKlein){
    fixMatrix = invGenerators[1];
    return true;
  }
  if(kleinSamplePoint.y > halfCubeWidthKlein){
    fixMatrix = invGenerators[2];
    return true;
  }
  if(kleinSamplePoint.y < -halfCubeWidthKlein){
    fixMatrix = invGenerators[3];
    return true;
  }
  if(kleinSamplePoint.z > halfCubeWidthKlein){
    fixMatrix = invGenerators[4];
    return true;
  }
  if(kleinSamplePoint.z < -halfCubeWidthKlein){
    fixMatrix = invGenerators[5];
    return true;
  }
  return false;
}

void raymarch(vec4 rO, vec4 rD){
  int fakeI = 0;
  float globalDepth = MIN_DIST;
  float localDepth = globalDepth;
  mat4 fixMatrix = mat4(1.0);
  vec4 localrO = rO;
  vec4 localrD = rD;
  
  // Trace the local scene, then the global scene:
  for(int i = 0; i< MAX_MARCHING_STEPS; i++){
    if(fakeI >= maxSteps || globalDepth >= maxDist){
      //when we break it's as if we reached our max marching steps
      break;
    }
    fakeI++;
    vec4 localEndPoint = pointOnGeodesic(localrO, localrD, localDepth);
    if(isOutsideCell(localEndPoint, fixMatrix)){
      totalFixMatrix *= fixMatrix;
      localrO = eucNormalize(localEndPoint*fixMatrix, false);
      localrD = eucFixDirection(localrO, localrD, fixMatrix); 
      localDepth = MIN_DIST;
    }
    else{
      float localDist = min(0.5,localSceneSDF(localEndPoint));
      if(localDist < EPSILON){
        hitWhich = 3;
        sampleEndPoint = localEndPoint;
        sampleTangentVector = tangentVectorOnGeodesic(localrO, localrD, localDepth);
        break;
      }
      localDepth += localDist;
      globalDepth += localDist;
    }
  }
  
  // Set localDepth to our new max tracing distance:
  localDepth = min(globalDepth, maxDist);
  globalDepth = MIN_DIST;
  fakeI = 0;
  for(int i = 0; i< MAX_MARCHING_STEPS; i++){
    if(fakeI >= maxSteps){
      break;
    }
    fakeI++;
    vec4 globalEndPoint = pointOnGeodesic(rO, rD, globalDepth);
    float globalDist = globalSceneSDF(globalEndPoint);
    if(globalDist < EPSILON){
      // hitWhich has been set by globalSceneSDF
      sampleEndPoint = globalEndPoint;
      sampleTangentVector = tangentVectorOnGeodesic(rO, rD, globalDepth);
      return;
    }
    globalDepth += globalDist;
    if(globalDepth >= localDepth){
      break;
    }
  }
}

void main(){
    vec4 rayOrigin = ORIGIN;

    //stereo translations ----------------------------------------------------
    bool isLeft = gl_FragCoord.x/screenResolution.x <= 0.5;
    vec4 rayDirV = getRayPoint(screenResolution, gl_FragCoord.xy, isLeft);
    if(isStereo == 1){
        if(isLeft){
            rayOrigin *= stereoBoosts[0];
            rayDirV *= stereoBoosts[0];
        }
        else{
            rayOrigin *= stereoBoosts[1];
            rayDirV *= stereoBoosts[1];
        }
    }

    //camera position must be translated in hyperboloid -----------------------
    rayOrigin *= currentBoost;
    rayDirV *= currentBoost;

    //generate direction then transform to hyperboloid ------------------------
    vec4 rayDirVPrime = eucDirection(rayOrigin, rayDirV);

    //get our raymarched distance back ------------------------
    raymarch(rayOrigin, rayDirVPrime);

    //Based on hitWhich decide whether we hit a global object, local object, or nothing
    if(hitWhich == 0){ //Didn't hit anything ------------------------
        gl_FragColor = vec4(0.0);
        return;
    }
    else if(hitWhich == 1){ // global lights
        gl_FragColor = vec4(globalLightColor.rgb, 1.0);
        return;
    }
    else{ // objects
        N = estimateNormal(sampleEndPoint);
        vec3 color;
        if(hitWhich == 2){ // global objects
            color = phongModel(invGlobalObjectBoosts[0], true);
        }else{ // local objects
            color = phongModel(mat4(1.0), false);
        }
        gl_FragColor = vec4(color, 1.0);
    }
}
END FRAGMENT