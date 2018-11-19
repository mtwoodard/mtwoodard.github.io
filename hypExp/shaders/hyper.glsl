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
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;
const vec4 ORIGIN = vec4(0,0,0,1);
//--------------------------------------------
//Generated Constants
//--------------------------------------------
const float halfIdealCubeWidthKlein = 0.5773502692;
const vec4 idealCubeCornerKlein = vec4(halfIdealCubeWidthKlein, halfIdealCubeWidthKlein, halfIdealCubeWidthKlein, 1.0);
//--------------------------------------------
//Global Variables
//--------------------------------------------
vec4 sampleEndPoint = vec4(1, 1, 1, 1);
vec4 sampleTangentVector = vec4(1, 1, 1, 1);
vec4 N = ORIGIN; //normal vector
vec4 globalLightColor = ORIGIN;
int hitWhich = 0;
//-------------------------------------------
//Translation & Utility Variables
//--------------------------------------------
uniform int isStereo;
uniform vec2 screenResolution;
uniform float fov;
uniform mat4 invGenerators[6];
uniform mat4 currentBoost;
uniform mat4 stereoBoosts[2];
uniform mat4 cellBoost; 
uniform mat4 invCellBoost;
uniform int maxSteps;
//--------------------------------------------
//Lighting Variables & Global Object Variables
//--------------------------------------------
uniform vec4 lightPositions[4];
//w component is the light's attenuation -- 6 since we need controllers
uniform vec4 lightIntensities[6]; 
uniform int attnModel;
uniform sampler2D texture;
//Max is two
//uniform int controllerCount; 
//uniform mat4 controllerBoosts[2];
//uniform vec4 controllerDualPoints[6];
uniform mat4 globalObjectBoosts[4];
uniform mat4 invGlobalObjectBoosts[4]; 
uniform vec3 globalObjectRadii[4];
uniform int globalObjectTypes[4];
//--------------------------------------------
//Scene Dependent Variables
//--------------------------------------------
uniform float halfCubeWidthKlein;
uniform float sphereRad;
uniform float horosphereSize;

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

float sinh(float x){
  float eX = exp(x);
  return (0.5 * (eX - 1.0/eX));
}

float asinh(float x){
  return log(x + sqrt(x*x+1.0));
}

//-------------------------------------------------------
// Hyperbolic Functions
//-------------------------------------------------------
float hypDot(vec4 u, vec4 v){
  return u.x*v.x + u.y*v.y + u.z*v.z - u.w*v.w; // hyp Dot
}
float hypNorm(vec4 v){
  return sqrt(abs(hypDot(v,v)));
}
vec4 hypNormalize(vec4 u, bool toTangent){
  return u/hypNorm(u);
}
float hypDistance(vec4 u, vec4 v){
  float bUV = -hypDot(u,v);
  return acosh(bUV);
}

//Given two positions find the unit tangent vector at the first that points to the second
vec4 hypDirection(vec4 u, vec4 v){
  vec4 w = v + hypDot(u,v)*u;
  return hypNormalize(w, true);
}

//calculate the new direction vector (v) for the continuation of the ray from the new ray origin (u)
//having moved by fix matrix
vec4 hypFixDirection(vec4 u, vec4 v, mat4 fixMatrix){
  return hypDirection(u, v*fixMatrix); 
}

//-------------------------------------------------------
//Hyperboloid Functions
//-------------------------------------------------------

vec4 projectToKlein(vec4 v){
  return v/v.w;
}

// Get point at distance dist on the geodesic from u in the direction vPrime
vec4 pointOnGeodesic(vec4 u, vec4 vPrime, float dist){
  return u*cosh(dist) + vPrime*sinh(dist);
}

vec4 tangentVectorOnGeodesic(vec4 u, vec4 vPrime, float dist){
  // note that this point has hypDot with itself of -1, so it is on other hyperboloid
  return u*sinh(dist) + vPrime*cosh(dist);
}

vec4 pointOnGeodesicAtInfinity(vec4 u, vec4 vPrime){ // returns point on the light
  // cone intersect Klein model corresponding to the point at infinity on the
  // geodesic through u and v
  return projectToKlein(u + vPrime);
}

//---------------------------------------------------------------------
// Sign Distance Fields
//---------------------------------------------------------------------

float sphereSDF(vec4 samplePoint, vec4 center, float radius){
  return hypDistance(samplePoint, center) - radius;
}

float sortOfEllipsoidSDF(vec4 samplePoint, mat4 boostMatrix){
  return sphereSDF(hypNormalize(samplePoint * boostMatrix, false), ORIGIN, 0.05);
}

// A horosphere can be constructed by offseting from a standard horosphere.
// Our standard horosphere will have a center in the direction of lightPoint
// and go through the origin. Negative offsets will 'shrink' it.
float horosphereHSDF(vec4 samplePoint, vec4 lightPoint, float offset){
  return log(-hypDot(samplePoint, lightPoint)) - offset;
}

float localSceneSDF(vec4 samplePoint){
    float sphere = sphereSDF(samplePoint, ORIGIN, sphereRad);
    float vertexSphere = horosphereHSDF(abs(samplePoint), idealCubeCornerKlein, horosphereSize);
    float final = -unionSDF(vertexSphere,sphere);
    return final;
}

//GLOBAL OBJECTS SCENE ++++++++++++++++++++++++++++++++++++++++++++++++
float globalSceneSDF(vec4 samplePoint){
  vec4 absoluteSamplePoint = samplePoint * cellBoost; // correct for the fact that we have been moving
  float distance = MAX_DIST;
  //Light Objects
  for(int i=0; i<NUM_LIGHTS; i++){
    float objDist;
    if(lightIntensities[i].w == 0.0) { objDist = MAX_DIST; }
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
  /*Controller Objects
  for(int i=0; i<2; i++){
    if(controllerCount != 0){
      //float objDist = sphereSDF(absoluteSamplePoint, ORIGIN*controllerBoosts[i-4]*currentBoost, 1.0/(10.0 * lightIntensities[i].w));
      float objDist = controllerSDF(absoluteSamplePoint, controllerBoosts[i-4]*currentBoost, 1.0/(10.0 * lightIntensities[i].w));
      distance = min(distance, objDist);
      if(distance < EPSILON){
        hitWhich = 1;
        globalLightColor = lightIntensities[i+4];
        return distance;
      }
      if(controllerCount == 1) break;
    }
  }*/
  //Global Objects
  for(int i=0; i<NUM_OBJECTS; i++) {
    float objDist;
    if(length(globalObjectRadii[i]) == 0.0){ objDist = MAX_DIST;}
    else{
      if(globalObjectTypes[i] == 0) { objDist = sphereSDF(absoluteSamplePoint, globalObjectBoosts[i][3], globalObjectRadii[i].x); }
      else if(globalObjectTypes[i] == 1) { objDist = sortOfEllipsoidSDF(absoluteSamplePoint, globalObjectBoosts[i]);}
      else { objDist = MAX_DIST; }
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
    vec3 n = hypNormalize(N*toOrigin, true).xyz; //Very hacky you are warned
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
  float distToLight = hypDistance(SP, TLP);
  float att = attenuation(distToLight, lightIntensity);
  //Calculations - Phong Reflection Model
  float shadow = 1.0;
  vec4 L = hypDirection(SP, TLP);
  vec4 R = 2.0*hypDot(L, N)*N - L;
  //Calculate Diffuse Component
  float nDotL = max(hypDot(N, L),0.0);
  vec3 diffuse = lightIntensity.rgb * nDotL;
  //Calculate Specular Component
  float rDotV = max(hypDot(R, V),0.0);
  vec3 specular = lightIntensity.rgb * pow(rDotV,10.0);
  //Compute final color
  return att*(shadow*((diffuse*baseColor) + specular));
}

vec3 phongModel(mat4 invObjectBoost, bool isGlobal, mat4 totalFixMatrix){
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
   /* for(int i = 0; i<2; i++){
      if(controllerCount == 0) break; //if there are no controllers do nothing
      else translatedLightPosition = ORIGIN*controllerBoosts[i]*currentBoost;
      color += lightingCalculations(samplePoint, translatedLightPosition, V, baseColor, lightIntensities[i+4]);
      if(controllerCount == 1) break; //if there is one controller only do one loop
    }*/
    return color;
}

//NORMAL FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++
vec4 estimateNormal(vec4 p) { // normal vector is in tangent hyperplane to hyperboloid at p
    // float denom = sqrt(1.0 + p.x*p.x + p.y*p.y + p.z*p.z);  // first, find basis for that tangent hyperplane
    float newEp = EPSILON * 10.0;
    vec4 basis_x = hypNormalize(vec4(p.w,0.0,0.0,p.x), true);  // dw/dx = x/w on hyperboloid
    vec4 basis_y = vec4(0.0,p.w,0.0,p.y);  // dw/dy = y/denom
    vec4 basis_z = vec4(0.0,0.0,p.w,p.z);  // dw/dz = z/denom  /// note that these are not orthonormal!
    basis_y = hypNormalize(basis_y - hypDot(basis_y, basis_x)*basis_x, true); // need to Gram Schmidt
    basis_z = hypNormalize(basis_z - hypDot(basis_z, basis_x)*basis_x - hypDot(basis_z, basis_y)*basis_y, true);
    if(hitWhich != 3){ //global light scene
        return hypNormalize( //p+EPSILON*basis_x should be lorentz normalized however it is close enough to be good enough
            basis_x * (globalSceneSDF(p + newEp*basis_x) - globalSceneSDF(p - newEp*basis_x)) +
            basis_y * (globalSceneSDF(p + newEp*basis_y) - globalSceneSDF(p - newEp*basis_y)) +
            basis_z * (globalSceneSDF(p + newEp*basis_z) - globalSceneSDF(p - newEp*basis_z)),
            true
        );
    }
    else{ //local scene
        return hypNormalize(
            basis_x * (localSceneSDF(p + newEp*basis_x) - localSceneSDF(p - newEp*basis_x)) +
            basis_y * (localSceneSDF(p + newEp*basis_y) - localSceneSDF(p - newEp*basis_y)) +
            basis_z * (localSceneSDF(p + newEp*basis_z) - localSceneSDF(p - newEp*basis_z)),
            true
        );
    }
}

vec4 getRayPoint(vec2 resolution, vec2 fragCoord, bool isLeft){ //creates a point that our ray will go through
    if(isStereo == 1){
      resolution.x = resolution.x * 0.5;
      if(!isLeft) { fragCoord.x = fragCoord.x - resolution.x; }
    }
    vec2 xy = 0.2*((fragCoord - 0.5*resolution)/resolution.x);
    float z = 0.1/tan(radians(fov*0.5));
    vec4 p =  hypNormalize(vec4(xy,-z,1.0), false);
    return vec4(xy,-z,1.0);
}

// This function is intended to be geometry-agnostic.
// We should update some of the variable names.
bool isOutsideCell(vec4 samplePoint, out mat4 fixMatrix){
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

void raymarch(vec4 rO, vec4 rD, out mat4 totalFixMatrix){
    int fakeI = 0;
    float globalDepth = MIN_DIST; float localDepth = globalDepth;
    vec4 localrO = rO; vec4 localrD = rD;
    mat4 fixMatrix = mat4(1.0);
    totalFixMatrix = mat4(1.0); 
  
    // Trace the local scene, then the global scene:
    for(int i = 0; i< MAX_MARCHING_STEPS; i++){
        if(fakeI >= maxSteps || globalDepth >= MAX_DIST){
            //when we break it's as if we reached our max marching steps
            break;
        }
        fakeI++;
        vec4 localEndPoint = pointOnGeodesic(localrO, localrD, localDepth);
        if(isOutsideCell(localEndPoint, fixMatrix)){
            totalFixMatrix *= fixMatrix;
            localrO = hypNormalize(localEndPoint*fixMatrix, false);
            localrD = hypFixDirection(localrO, localrD, fixMatrix); 
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
    localDepth = min(globalDepth, MAX_DIST);
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
            // hitWhich has now been set
            totalFixMatrix = mat4(1.0);
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
    /*if(isStereo == 1){
        if(isLeft){
            rayOrigin *= stereoBoosts[0];
            rayDirV *= stereoBoosts[0];
        }
        else{
            rayOrigin *= stereoBoosts[1];
            rayDirV *= stereoBoosts[1];
        }
    }*/

    //camera position must be translated in hyperboloid -----------------------
    rayOrigin *= currentBoost;
    rayDirV *= currentBoost;
    //generate direction then transform to hyperboloid ------------------------
    //vec4 rayDirVPrime = hypDirection(rayOrigin, rayDirV);
    //get our raymarched distance back ------------------------
   // mat4 totalFixMatrix = mat4(1.0);
    //raymarch(rayOrigin, rayDirVPrime, totalFixMatrix);
    //gl_FragColor = vec4(0.5,0.0,0.0,1.0);
    //Based on hitWhich decide whether we hit a global object, local object, or nothing
    gl_FragColor = rayDirV;
    /*if(hitWhich == 0){ //Didn't hit anything ------------------------
        gl_FragColor = vec4(0.1);
        return;
    }
    else if(hitWhich == 1){ // global lights
        gl_FragColor = vec4(globalLightColor.rgb, 1.0);
        //gl_FragColor = vec4(0.0, 1.0, 0.0, 1.0);
        return;
    }
    else{ // objects
      N = estimateNormal(sampleEndPoint);
      vec3 color;
      if(hitWhich == 2){ // global objects
        color = phongModel(invGlobalObjectBoosts[0], true, totalFixMatrix);
      }else{ // local objects
        color = phongModel(mat4(1.0), false, totalFixMatrix);
      }
      gl_FragColor = vec4(color, 1.0);
    }*/
}
END FRAGMENT