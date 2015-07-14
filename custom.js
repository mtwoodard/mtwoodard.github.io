var app = angular.module('mtwoodardApp', []);
app.controller('myCtrl', function($scope){
        $scope.secret = false;
        $scope.pumpedUp = false;
        var koCode = [38, 38, 40, 40, 37, 39, 37, 39, 66, 65, 13];
        var count = 0;
        document.onkeydown = function(event){
          var key_code = event.keyCode;
          if(key_code === koCode[count]){
              count++;
              if(count === koCode.length){
                $scope.secret = true;
                console.log($scope.secret);
              }
          }
          else{
              count = 0;
          }
        }

        $scope.pumpIt = function(){
          $scope.pumpedUp = true;
          console.log('ayyyyyy');
        };
});
