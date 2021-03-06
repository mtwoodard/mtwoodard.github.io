BEGIN VERTEX
  void main()
  {
      gl_Position = projectionMatrix * modelViewMatrix * vec4(position.xyz, 1.0);
  }
END VERTEX

BEGIN FRAGMENT
  //--------------------------------------------
  //Global Constants
  //--------------------------------------------
  const int MAX_MARCHING_STEPS = 32;
  const float MIN_DIST = 0.0;
  const float MAX_DIST = 100.0;
  const float EPSILON = 0.0001;
  const float fov = 90.0;
  const float horosphereSize = -0.951621;
  const float sphereRad = 0.996216;
  const float halfCubeWidthKlein = 0.5773502692;
  const float globalObjectRadius = 0.2;
  const vec4 ORIGIN = vec4(0,0,0,1);
  //--------------------------------------------
  //Generated Constants
  //--------------------------------------------
  const float halfIdealCubeWidthKlein = 0.5773502692;
  const vec4 idealCubeCornerKlein = vec4(halfIdealCubeWidthKlein, halfIdealCubeWidthKlein, halfIdealCubeWidthKlein, 1.0);
  //-------------------------------------------
  //Translation & Utility Variables
  //--------------------------------------------
  uniform vec2 screenResolution;
  uniform mat4 invGenerators[6];
  uniform mat4 currentBoost;
  uniform mat4 cellBoost; 
  uniform mat4 invCellBoost;
  //--------------------------------------------
  //Lighting Variables & Global Object Variables
  //--------------------------------------------
  uniform vec4 lightPositions[4];
  uniform vec4 lightIntensities[4];
  
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
  
  float hypDot(vec4 u, vec4 v){
    return u.x*v.x + u.y*v.y + u.z*v.z - u.w*v.w; // Lorentz Dot
  }
  
  float hypNorm(vec4 v){
    return sqrt(abs(hypDot(v,v)));
  }
  
  vec4 hypNormalize(vec4 u){
    return u/hypNorm(u);
  }
  
  float hypDistance(vec4 u, vec4 v){
    float bUV = -hypDot(u,v);
    return acosh(bUV);
  }
  
  vec4 hypDirection(vec4 u, vec4 v){
    vec4 w = v + hypDot(u,v)*u;
    return hypNormalize(w);
  }
  
  //-------------------------------------------------------
  //Hyperboloid Functions
  //-------------------------------------------------------
  // Get point at distance dist on the geodesic from u in the direction vPrime
  vec4 pointOnGeodesic(vec4 u, vec4 vPrime, float dist){
    return u*cosh(dist) + vPrime*sinh(dist);
  }
  
  vec4 tangentVectorOnGeodesic(vec4 u, vec4 vPrime, float dist){
    // note that this point has hypDot with itself of -1, so it is on other hyperboloid
    return u*sinh(dist) + vPrime*cosh(dist);
  }
  
  //---------------------------------------------------------------------
  //Raymarch Primitives
  //---------------------------------------------------------------------
  // A horosphere can be constructed by offseting from a standard horosphere.
  // Our standard horosphere will have a center in the direction of lightPoint
  // and go through the origin. Negative offsets will shrink it.
  float horosphereHSDF(vec4 samplePoint, vec4 lightPoint, float offset){
    return log(-hypDot(samplePoint, lightPoint)) - offset;
  }
  
  float sphereSDF(vec4 samplePoint, vec4 center, float radius){
    return hypDistance(samplePoint, center) - radius;
  }
  
  //---------------------------------------------------------------------
  //Scene Definitions
  //---------------------------------------------------------------------
  float localSceneSDF(vec4 samplePoint){
    float sphere = sphereSDF(samplePoint, ORIGIN, sphereRad);
    float vertexSphere = 0.0;
    vertexSphere = horosphereHSDF(abs(samplePoint), idealCubeCornerKlein, horosphereSize);
    float final = -min(vertexSphere,sphere); //unionSDF
    return final;
  }
  
  //--------------------------------------------------------------------
  // Lighting Functions
  //--------------------------------------------------------------------
  vec3 lightingCalculations(vec4 SP, vec4 TLP, vec4 V, vec4 normal, vec4 lightIntensity){
    //Calculations - Phong Reflection Model
    vec4 L = hypDirection(SP, TLP);
    vec4 R = 2.0*hypDot(L, normal)*normal - L;
    //Calculate Diffuse Component
    float nDotL = max(hypDot(normal, L),0.0);
    vec3 diffuse = lightIntensity.rgb * nDotL;
    //Calculate Specular Component
    float rDotV = max(hypDot(R, V),0.0);
    vec3 specular = lightIntensity.rgb * pow(rDotV,10.0);
    //Attenuation - Inverse Square
    float distToLight = hypDistance(SP, TLP);
    float att = 1.0/(0.01 + lightIntensity.w * distToLight* distToLight);
    //Compute final color
    return att*(diffuse + specular);
  }
  
  vec3 phongModel(vec4 samplePoint, vec4 tangentVector, vec4 normal, mat4 totalFixMatrix){
    vec4 V = -tangentVector;
    vec3 color = vec3(0.0);
    //--------------------------------------------
    //Lighting Calculations
    //--------------------------------------------
    vec4 translatedLightPosition = vec4(0.0);
    //Standard Light Objects
    for(int i = 0; i<4; i++){ //4 is the number of lights we can use
      if(lightIntensities[i].w != 0.0){
        translatedLightPosition = lightPositions[i]*invCellBoost*totalFixMatrix;
        color += lightingCalculations(samplePoint, translatedLightPosition, V, normal, lightIntensities[i]);
      }
    }
    return color;
  }
  
  
  //NORMAL FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++
  vec4 estimateNormal(vec4 p) { // normal vector is in tangent hyperplane to hyperboloid at p
    float newEp = EPSILON * 10.0;
    vec4 basis_x = hypNormalize(vec4(p.w,0.0,0.0,p.x));  // dw/dx = x/w on hyperboloid
    vec4 basis_y = vec4(0.0,p.w,0.0,p.y);  // dw/dy = y/denom
    vec4 basis_z = vec4(0.0,0.0,p.w,p.z);  // dw/dz = z/denom  /// note that these are not orthonormal!
    basis_y = hypNormalize(basis_y - hypDot(basis_y, basis_x)*basis_x); // need to Gram Schmidt
    basis_z = hypNormalize(basis_z - hypDot(basis_z, basis_x)*basis_x - hypDot(basis_z, basis_y)*basis_y);
    return hypNormalize(
       basis_x * (localSceneSDF(p + newEp*basis_x) - localSceneSDF(p - newEp*basis_x)) +
      basis_y * (localSceneSDF(p + newEp*basis_y) - localSceneSDF(p - newEp*basis_y)) +
      basis_z * (localSceneSDF(p + newEp*basis_z) - localSceneSDF(p - newEp*basis_z)));
  }
  
  vec4 getRayPoint(vec2 resolution, vec2 fragCoord){ //creates a point that our ray will go through
    vec2 xy = 0.2*((fragCoord - 0.5*resolution)/resolution.x);
    float z = 0.1/tan(radians(fov*0.5));
    vec4 p =  hypNormalize(vec4(xy,-z,1.0));
    return p;
  }

  bool isOutsideCell(vec4 samplePoint, out mat4 fixMatrix){
    vec4 kleinSamplePoint = samplePoint/samplePoint.w; //project to klein
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
  
  bool raymarch(vec4 rO, vec4 rD, out vec4 samplePoint, out vec4 tangent, out mat4 totalFixMatrix){
    float globalDepth = MIN_DIST; float localDepth = globalDepth;
    totalFixMatrix = mat4(1.0);
    for(int i = 0; i< MAX_MARCHING_STEPS; i++){
      mat4 fixMatrix; //= mat4(1.0);
      vec4 endPoint = pointOnGeodesic(rO, rD, localDepth);
      if(isOutsideCell(endPoint, fixMatrix)){
        totalFixMatrix *= fixMatrix;
        rO = hypNormalize(endPoint*fixMatrix);
        rD = hypDirection(rO, rD*fixMatrix);
        localDepth = MIN_DIST;
      }
      else{
        float dist = localSceneSDF(endPoint); 
        if(dist < EPSILON){
          //Pass information out to global variables
          samplePoint = endPoint; //local sample point
          tangent = tangentVectorOnGeodesic(rO, rD, localDepth); //local tangent vector
          return true;
        }
        globalDepth += dist;
        localDepth += dist;
        if(globalDepth >= MAX_DIST) return false;
      }
    }
    return false;
  }
  
  void main(){
    //get intial origin and ray direction -------------------------------------
    vec4 rayOrigin = ORIGIN;
    vec4 rayDirV = getRayPoint(screenResolution, gl_FragCoord.xy);
    //camera position must be translated in hyperboloid -----------------------
    rayOrigin *= currentBoost;
    rayDirV *= currentBoost;
    //generate direction then transform to hyperboloid ------------------------
     vec4 rayDirVPrime = hypDirection(rayOrigin, rayDirV);
    //get our raymarch info ---------------------------------------------------
    mat4 totalFixMatrix = mat4(1.0);
    vec4 samplePoint; vec4 tangent;
    bool hit = raymarch(rayOrigin, rayDirVPrime, samplePoint, tangent, totalFixMatrix);
    //if we hit something color it --------------------------------------------
    if(!hit) gl_FragColor = vec4(0.05,0.05,0.05, 1.0);
    else{
      vec4 normal = estimateNormal(samplePoint);
      vec3 color = phongModel(samplePoint, tangent, normal, totalFixMatrix);
      gl_FragColor = vec4(color, 1.0);
    }
  }
END FRAGMENT