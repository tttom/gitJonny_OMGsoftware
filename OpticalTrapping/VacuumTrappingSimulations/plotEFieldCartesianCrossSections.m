function plotEFieldCartesianCrossSections(xRange,yRange,zRange,psf,psfField)
    
    % Determine closest-to-zero pixel coords
    [x_min_val x_min_pos] = min(abs(xRange));
    [y_min_val y_min_pos] = min(abs(yRange));
    [z_min_val z_min_pos] = min(abs(zRange));
    
    % Take 2D cross-sections
    psf_xy = psf(:,:,z_min_pos);
    psf_xz = psf(:,y_min_pos,:);
    psf_yz = psf(x_min_pos,:,:);
    
    psfField_xy = psfField(:,:,z_min_pos,:);
    psfField_xz = psfField(:,y_min_pos,:,:);
    psfField_yz = psfField(x_min_pos,:,:,:);
    
    % Determine max field amplitude
    maxFieldAmp = max(max(max(abs(psfField_xy(:))),max(abs(psfField_xz(:)))),max(abs(psfField_yz(:))));
    
    figure;
    subplot(4,6,[1 2]);imagesc(xRange * 1e6,yRange * 1e6,psf_xy.');axis image;
    xlabel('x-axis [um]');ylabel('y-axis [um]');
    title('Intensity (z=0)');
    subplot(4,6,7);imagesc(xRange * 1e6,yRange * 1e6,abs(squeeze(psfField_xy(:,:,1,1))).');axis image;
    xlabel('x-axis [um]');ylabel('y-axis [um]');
    title('abs(E_x)');
    set(gca,'CLim',[0 maxFieldAmp]);
    subplot(4,6,8);imagesc(xRange * 1e6,yRange * 1e6,angle(squeeze(psfField_xy(:,:,1,1))).');axis image;
    xlabel('x-axis [um]');ylabel('y-axis [um]');
    title('angle(E_x)');
    subplot(4,6,13);imagesc(xRange * 1e6,yRange * 1e6,abs(squeeze(psfField_xy(:,:,1,2))).');axis image;
    xlabel('x-axis [um]');ylabel('y-axis [um]');
    title('abs(E_y)');
    set(gca,'CLim',[0 maxFieldAmp]);
    subplot(4,6,14);imagesc(xRange * 1e6,yRange * 1e6,angle(squeeze(psfField_xy(:,:,1,2))).');axis image;
    xlabel('x-axis [um]');ylabel('y-axis [um]');
    title('angle(E_y)');
    subplot(4,6,19);imagesc(xRange * 1e6,yRange * 1e6,abs(squeeze(psfField_xy(:,:,1,3))).');axis image;
    xlabel('x-axis [um]');ylabel('y-axis [um]');
    title('abs(E_z)');
    set(gca,'CLim',[0 maxFieldAmp]);
    subplot(4,6,20);imagesc(xRange * 1e6,yRange * 1e6,angle(squeeze(psfField_xy(:,:,1,3))).');axis image;
    xlabel('x-axis [um]');ylabel('y-axis [um]');
    title('angle(E_z)');
    
    subplot(4,6,[3 4]);imagesc(xRange * 1e6,zRange * 1e6,squeeze(psf_xz).');axis image;
    xlabel('x-axis [um]');ylabel('z-axis [um]');
    title('Intensity (y=0)');
    subplot(4,6,9);imagesc(xRange * 1e6,zRange * 1e6,abs(squeeze(psfField_xz(:,1,:,1))).');axis image;
    xlabel('x-axis [um]');ylabel('z-axis [um]');
    title('abs(E_x)');
    set(gca,'CLim',[0 maxFieldAmp]);
    subplot(4,6,10);imagesc(xRange * 1e6,zRange * 1e6,angle(squeeze(psfField_xz(:,1,:,1))).');axis image;
    xlabel('x-axis [um]');ylabel('z-axis [um]');
    title('angle(E_x)');
    subplot(4,6,15);imagesc(xRange * 1e6,zRange * 1e6,abs(squeeze(psfField_xz(:,1,:,2))).');axis image;
    xlabel('x-axis [um]');ylabel('z-axis [um]');
    title('abs(E_y)');
    set(gca,'CLim',[0 maxFieldAmp]);
    subplot(4,6,16);imagesc(xRange * 1e6,zRange * 1e6,angle(squeeze(psfField_xz(:,1,:,2))).');axis image;
    xlabel('x-axis [um]');ylabel('z-axis [um]');
    title('angle(E_y)');
    subplot(4,6,21);imagesc(xRange * 1e6,zRange * 1e6,abs(squeeze(psfField_xz(:,1,:,3))).');axis image;
    xlabel('x-axis [um]');ylabel('z-axis [um]');
    title('abs(E_z)');
    set(gca,'CLim',[0 maxFieldAmp]);
    subplot(4,6,22);imagesc(xRange * 1e6,zRange * 1e6,angle(squeeze(psfField_xz(:,1,:,3))).');axis image;
    xlabel('x-axis [um]');ylabel('z-axis [um]');
    title('angle(E_z)');
    
    subplot(4,6,[5 6]);imagesc(yRange * 1e6,zRange * 1e6,squeeze(psf_yz).');axis image;
    xlabel('y-axis [um]');ylabel('z-axis [um]');
    title('Intensity (x=0)');
    subplot(4,6,11);imagesc(yRange * 1e6,zRange * 1e6,abs(squeeze(psfField_yz(1,:,:,1))).');axis image;
    xlabel('y-axis [um]');ylabel('z-axis [um]');
    title('abs(E_x)');
    set(gca,'CLim',[0 maxFieldAmp]);
    subplot(4,6,12);imagesc(yRange * 1e6,zRange * 1e6,angle(squeeze(psfField_yz(1,:,:,1))).');axis image;
    xlabel('y-axis [um]');ylabel('z-axis [um]');
    title('angle(E_x)');
    subplot(4,6,17);imagesc(yRange * 1e6,zRange * 1e6,abs(squeeze(psfField_yz(1,:,:,2))).');axis image;
    xlabel('y-axis [um]');ylabel('z-axis [um]');
    title('abs(E_y)');
    set(gca,'CLim',[0 maxFieldAmp]);
    subplot(4,6,18);imagesc(yRange * 1e6,zRange * 1e6,angle(squeeze(psfField_yz(1,:,:,2))).');axis image;
    xlabel('y-axis [um]');ylabel('z-axis [um]');
    title('angle(E_y)');
    subplot(4,6,23);imagesc(yRange * 1e6,zRange * 1e6,abs(squeeze(psfField_yz(1,:,:,3))).');axis image;
    xlabel('y-axis [um]');ylabel('z-axis [um]');
    title('abs(E_z)');
    set(gca,'CLim',[0 maxFieldAmp]);
    subplot(4,6,24);imagesc(yRange * 1e6,zRange * 1e6,angle(squeeze(psfField_yz(1,:,:,3))).');axis image;
    xlabel('y-axis [um]');ylabel('z-axis [um]');
    title('angle(E_z)');
    
end