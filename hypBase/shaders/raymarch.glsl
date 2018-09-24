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
  const lowp int MAX_MARCHING_STEPS = 32;
  const lowp float MIN_DIST = 0.0;
  const lowp float MAX_DIST = 100.0;
  const lowp float EPSILON = 0.0001;
  const lowp float fov = 90.0;
  const lowp float horosphereSize = -0.951621;
  const lowp float sphereRad = 0.996216;
  const lowp float halfCubeWidthKlein = 0.5773502692;
  const lowp float globalObjectRadius = 0.2;
  const lowp vec4 ORIGIN = vec4(0,0,0,1);
  //--------------------------------------------
  //Generated Constants
  //--------------------------------------------
  const lowp float halfIdealCubeWidthKlein = 0.5773502692;
  const lowp vec4 idealCubeCornerKlein = vec4(halfIdealCubeWidthKlein, halfIdealCubeWidthKlein, halfIdealCubeWidthKlein, 1.0);
  //-------------------------------------------
  //Translation & Utility Variables
  //--------------------------------------------
  uniform lowp vec2 screenResolution;
  uniform lowp mat4 invGenerators[6];
  uniform lowp mat4 currentBoost;
  uniform lowp mat4 cellBoost; 
  uniform lowp mat4 invCellBoost;
  //--------------------------------------------
  //Lighting Variables & Global Object Variables
  //--------------------------------------------
  uniform lowp vec4 lightPositions[4];
  uniform lowp vec4 lightIntensities[4];
  
  //--------------------------------------------------------------------
  // Hyperbolic Functions
  //--------------------------------------------------------------------
  lowp float cosh(lowp float x){
    lowp float eX = exp(x);
    return (0.5 * (eX + 1.0/eX));
  }
  
  lowp float acosh(lowp float x){ //must be more than 1
    return log(x + sqrt(x*x-1.0));
  }
  
  lowp float sinh(lowp float x){
    lowp float eX = exp(x);
    return (0.5 * (eX - 1.0/eX));
  }
  
  lowp float hypDot(lowp vec4 u, lowp vec4 v){
    return u.x*v.x + u.y*v.y + u.z*v.z - u.w*v.w; // Lorentz Dot
  }
  
  lowp float hypNorm(lowp vec4 v){
    return sqrt(abs(hypDot(v,v)));
  }
  
  lowp vec4 hypNormalize(lowp vec4 u){
    return u/hypNorm(u);
  }
  
  lowp float hypDistance(lowp vec4 u, lowp vec4 v){
    lowp float bUV = -hypDot(u,v);
    return acosh(bUV);
  }
  
  lowp vec4 hypDirection(lowp vec4 u, lowp vec4 v){
    lowp vec4 w = v + hypDot(u,v)*u;
    return hypNormalize(w);
  }
  
  //-------------------------------------------------------
  //Hyperboloid Functions
  //-------------------------------------------------------
  // Get point at distance dist on the geodesic from u in the direction vPrime
  lowp vec4 pointOnGeodesic(lowp vec4 u, lowp vec4 vPrime, lowp float dist){
    return u*cosh(dist) + vPrime*sinh(dist);
  }
  
  lowp vec4 tangentVectorOnGeodesic(lowp vec4 u, lowp vec4 vPrime, lowp float dist){
    // note that this point has hypDot with itself of -1, so it is on other hyperboloid
    return u*sinh(dist) + vPrime*cosh(dist);
  }
  
  //---------------------------------------------------------------------
  //Raymarch Primitives
  //---------------------------------------------------------------------
  // A horosphere can be constructed by offseting from a standard horosphere.
  // Our standard horosphere will have a center in the direction of lightPoint
  // and go through the origin. Negative offsets will shrink it.
  float horosphereHSDF(lowp vec4 samplePoint, lowp vec4 lightPoint, lowp float offset){
    return log(-hypDot(samplePoint, lightPoint)) - offset;
  }
  
  lowp float sphereSDF(lowp vec4 samplePoint, lowp vec4 center, lowp float radius){
    return hypDistance(samplePoint, center) - radius;
  }
  
  //---------------------------------------------------------------------
  //Scene Definitions
  //---------------------------------------------------------------------
  lowp float localSceneSDF(lowp vec4 samplePoint){
    lowp float sphere = sphereSDF(samplePoint, ORIGIN, sphereRad);
    lowp float vertexSphere = 0.0;
    vertexSphere = horosphereHSDF(abs(samplePoint), idealCubeCornerKlein, horosphereSize);
    lowp float final = -min(vertexSphere,sphere); //unionSDF
    return final;
  }
  
  //--------------------------------------------------------------------
  // Lighting Functions
  //--------------------------------------------------------------------
  lowp vec3 lightingCalculations(lowp vec4 SP, lowp vec4 TLP, lowp vec4 V, lowp vec4 normal, lowp vec4 lightIntensity){
    //Calculations - Phong Reflection Model
    lowp vec4 L = hypDirection(SP, TLP);
    lowp vec4 R = 2.0*hypDot(L, normal)*normal - L;
    //Calculate Diffuse Component
    lowp float nDotL = max(hypDot(normal, L),0.0);
    lowp vec3 diffuse = lightIntensity.rgb * nDotL;
    lowp //Calculate Specular Component
    float rDotV = max(hypDot(R, V),0.0);
    lowp vec3 specular = lightIntensity.rgb * pow(rDotV,10.0);
    //Attenuation - Inverse Square
    lowp float distToLight = hypDistance(SP, TLP);
    lowp float att = 1.0/(0.01 + lightIntensity.w * distToLight* distToLight);
    //Compute final color
    return att*((diffuse*vec3(1.0)) + specular);
  }
  
  lowp vec3 phongModel(lowp vec4 samplePoint, lowp vec4 tangentVector, lowp vec4 normal, lowp mat4 totalFixMatrix){
    lowp vec4 V = -tangentVector;
    lowp vec3 color = vec3(0.0);
    //--------------------------------------------
    //Lighting Calculations
    //--------------------------------------------
    vec4 translatedLightPosition;
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
  lowp vec4 estimateNormal(lowp vec4 p) { // normal vector is in tangent hyperplane to hyperboloid at p
    lowp float newEp = EPSILON * 10.0;
    lowp vec4 basis_x = hypNormalize(vec4(p.w,0.0,0.0,p.x));  // dw/dx = x/w on hyperboloid
    lowp vec4 basis_y = vec4(0.0,p.w,0.0,p.y);  // dw/dy = y/denom
    lowp vec4 basis_z = vec4(0.0,0.0,p.w,p.z);  // dw/dz = z/denom  /// note that these are not orthonormal!
    basis_y = hypNormalize(basis_y - hypDot(basis_y, basis_x)*basis_x); // need to Gram Schmidt
    basis_z = hypNormalize(basis_z - hypDot(basis_z, basis_x)*basis_x - hypDot(basis_z, basis_y)*basis_y);
    return hypNormalize(
       basis_x * (localSceneSDF(p + newEp*basis_x) - localSceneSDF(p - newEp*basis_x)) +
      basis_y * (localSceneSDF(p + newEp*basis_y) - localSceneSDF(p - newEp*basis_y)) +
      basis_z * (localSceneSDF(p + newEp*basis_z) - localSceneSDF(p - newEp*basis_z)));
  }
  
  lowp vec4 getRayPoint(lowp vec2 resolution, lowp vec2 fragCoord){ //creates a point that our ray will go through
    lowp vec2 xy = 0.2*((fragCoord - 0.5*resolution)/resolution.x);
    lowp float z = 0.1/tan(radians(fov*0.5));
    lowp vec4 p =  hypNormalize(vec4(xy,-z,1.0));
    return p;
  }

  bool isOutsideCell(lowp vec4 samplePoint, out lowp mat4 fixMatrix){
    lowp vec4 kleinSamplePoint = samplePoint/samplePoint.w; //project to klein
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
  
  bool raymarch(lowp vec4 rO, lowp vec4 rD, out lowp vec4 samplePoint, out lowp vec4 tangent, out lowp mat4 totalFixMatrix){
    lowp float globalDepth = MIN_DIST; lowp float localDepth = globalDepth;
    for(lowp int i = 0; i< MAX_MARCHING_STEPS; i++){
      lowp mat4 fixMatrix;
      lowp vec4 endPoint = pointOnGeodesic(rO, rD, localDepth);
      if(isOutsideCell(endPoint, fixMatrix)){
        totalFixMatrix *= fixMatrix;
        rO = hypNormalize(endPoint*fixMatrix);
        rD = hypDirection(rO, rD*fixMatrix);
        localDepth = MIN_DIST;
      }
      else{
        lowp float dist = localSceneSDF(endPoint); 
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
    lowp vec4 rayOrigin = ORIGIN;
    lowp vec4 rayDirV = getRayPoint(screenResolution, gl_FragCoord.xy);
    //camera position must be translated in hyperboloid -----------------------
    rayOrigin *= currentBoost;
    rayDirV *= currentBoost;
    //generate direction then transform to hyperboloid ------------------------
    lowp vec4 rayDirVPrime = hypDirection(rayOrigin, rayDirV);
    //get our raymarch info ---------------------------------------------------
    lowp mat4 totalFixMatrix = mat4(1.0);
    lowp vec4 samplePoint; lowp vec4 tangent;
    bool hit = raymarch(rayOrigin, rayDirVPrime, samplePoint, tangent, totalFixMatrix);
    //if we hit something color it --------------------------------------------
    if(!hit) gl_FragColor = vec4(0.0);
    else{
      lowp vec4 normal = estimateNormal(samplePoint);
      lowp vec3 color = phongModel(samplePoint, tangent, normal, totalFixMatrix);
      gl_FragColor = vec4(color, 1.0);
    }
  }
END FRAGMENT